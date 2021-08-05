"""
Kepler multiporcessing

Example:
from kepler.process import Proxy

"""

from types import MethodType

from .plot import kepplots
from .plot.manager import PlotManager

from .code import Kepler, KeplerData
from .datainterface import DataInterface

from .code.kepbin import _kepler

from multiprocessing import Process, Queue, JoinableQueue
import numpy as np

class RemoteError(Exception):
    pass
class CommandNotFoundError(RemoteError):
    def __str__(self):
        return ' [REMOTE] ERROR: commnand "{}" not found'.format(*self.args)
    __repr__ = __str__
class ExecutionError(RemoteError):
    def __str__(self):
        return ' [REMOTE] ERROR: executing {}(*{}, **{}). {}'.format(*self.args)
class ActionError(RemoteError):
    def __str__(self):
        return ' [REMOTE] ERROR: unknown action "{}".'.format(*self.args)
class MissingAttributeError(RemoteError):
    def __str__(self):
        return ' [REMOTE] ERROR: missing attribute'.format(*self.args)
class MissingValueError(RemoteError):
    def __str__(self):
        return ' [REMOTE] ERROR: missing value'.format(*self.args)
class RemoteAttributeError(RemoteError):
    def __str__(self):
        return ' [REMOTE] ERROR: unknown attribute "{}".'.format(*self.args)
class RemoteMethodError(RemoteError):
    def __str__(self):
        return ' [REMOTE] ERROR: trying to return remote method "{}".'.format(*self.args)
class RemoteIndexError(RemoteError):
    def __str__(self):
        return ' [REMOTE] ERROR: invalid index {} for attribute "{}".'.format(*self.args)
class RemoteException(RemoteError):
    def __str__(self):
        return ' [REMOTE] ERROR: exception {}.'.format(*self.args)

void = object()

from functools import partial

class KepProcess(Process):
    def __init__(self, task_queue, done_queue, call_queue, repl_queue, params):
        super().__init__()
        args = params.get('args', ())
        kwargs = params.get('kwargs', {})
        kwargs.setdefault('plot', False)
        self.kepler = Kepler(*args, **kwargs)
        self.task_queue = task_queue
        self.done_queue = done_queue
        self.call_queue = call_queue
        self.repl_queue = repl_queue

    def endkepler(self, code = None):
        print(' [DEBUG] remote end callback')
        task = {'action' : 'call',
                'func' : 'terminate',
                'kwargs' : {'code' : code}
                }
        self.call_queue.put(task)
        # repl = self.repl_queue.get()
        # print(f' [KepProcess] : {repl}')

    def plotkepler(self):
        print(' [DEBUG] remote plot callback')
        task = {'action' : 'call',
                'func' : 'plot',
                'kwargs' : {'itvstart' : self.kepler.kd.parm.itvstart,
                            'ipixtype' : self.kepler.kd.parm.ipixtype,
                            },
                }
        self.call_queue.put(task)
        # repl = self.repl_queue.get()
        # print(f' [KepProcess] : {repl}')

    def getskepler(self, msg):
        print(' [DEBUG] remote input callback')
        task = {'action' : 'call',
                'func' : 'gets',
                }
        self.call_queue.put(task)
        # repl = self.repl_queue.get()
        # print(f' [KepProcess] : {repl}')

        # s = repl['data']
        # msg[:] = 0
        # for i, c in enumerate(s):
        #     msg[i] = bytes(c, encoding = 'us-ascii')[0]

    def run(self):
        _kepler.getskepler = self.getskepler
        _kepler.endkepler = self.endkepler
        _kepler.plotkepler = self.plotkepler

        # this can be geeneralized easily
        for item in iter(self.task_queue.get, 'STOP'):
            self.done_queue.put(self.process(item))
            self.task_queue.task_done()
        self.task_queue.task_done()

    def __call__(self, params):
        return self.process(params)

    def process(self, params):
        action = params.get('action', 'call')
        if action == 'call':
            func = params.get('func', None)
            args = params.get('args', ())
            kwargs = params.get('kwargs', {})
            return self.call(func, args, kwargs)
        elif action in ('get', 'type'):
            attr = params.get('attr', None)
            if attr is None:
                return MissingAttributeError()
            attrs = attr.split('.')
            try:
                val = getattr(self.kepler, attrs[0])
            except AttributeError:
                return RemoteAttributeError(attr)
            for a in attrs[1:]:
                try:
                    val = getattr(val, a)
                except AttributeError:
                    return RemoteAttributeError(attr)
            if isinstance(val, MethodType):
                return RemoteMethodError(attr)
            if 'index' in params:
                index = params.get('index')
                try:
                    val = val[index]
                except:
                    return RemoteIndexError(index, attr)
            if action == 'get':
                return val
            else:
                return type(val)
        elif action == 'set':
            attr = params.get('attr', None)
            if attr is None:
                return MissingAttributeError()
            try:
                value = params.get('value')
            except:
                return MissingValueError()
            obj = self.kepler
            attrs = attr.split('.')
            for a in attrs[:-1]:
                try:
                    obj = getattr(obj, a)
                except AttributeError:
                    return RemoteAttributeError(attr)
            if 'index' in params:
                index = params.get('index')
                try:
                    var = getattr(obj, attrs[-1])
                except AttributeError:
                    return RemoteAttributeError(attr)
                try:
                    var[index] = value
                except IndexError:
                    return RemoteIndexError(index, attr)
                except Exception as error:
                    return RemoteException(error)
                return None
            else:
                try:
                    return setattr(obj, attrs[-1], value)
                except AttributeError:
                    return RemoteAttributeError(attr)
                except Exception as error:
                    return RemoteException(error)
        elif action == 'data':
            if not 'data' in self.__dict__:
                self.data = KeplerData()
            attr = params.get('attr', None)
            if attr is None:
                return MissingAttributeError()
            try:
                val = getattr(self.data, attr)
            except AttributeError:
                return RemoteAttributeError(attr)
            return val
        elif action == 'execute':
            cmd = params.get('cmd', None)
            if cmd is None:
                raise MissingCommandError()
            for c in cmd.splitlines():
                self.kepler.execute(c)
            return
        return ActionError(action)

    def call(self, func, args, kwargs):
        kwargs = dict(kwargs)
        kwargs['remote'] = True
        try:
            method = getattr(self.kepler, func)
        except:
            return CommandNotFoundError(func)
        try:
            result = method(*args, **kwargs)
        except Exception as error:
            return ExecutionError(func, args, kwargs, error)
        return result

class ProxyTemplate():
    _slots = (
        'remote',
        'task_queue',
        'done_queue',
        'call_queue',
        'repl_queue',
        'remote_arrays',
        'kp',
        'debug',
        'data',
        '_actions',
            )
    def __init__(self, *args, **kwargs):
        self.remote_arrays = kwargs.pop('remote_arrays', False)
        self.debug = kwargs.pop('debug', True)
        self._actions = list()
        params = dict(
            args = args,
            kwargs = kwargs,
            )
        self.setup_process(params)
        self.kp = PlotManager(self.data)

    def setup_process(self, args):
        """
        Needs to set up:

          self.remote
          self.task_queue
          self.done_queue
          self.call_queue
          self.repl_queue
        """
        raise NotImplementedError(f'Need to implement {__qualname__}')

    def __call__(self, func, *args, **kwargs):
        return self.call(func, *args, **kwargs)

    def handle_calls(self, action):
        pass
    # how to make sure there is no compiting calls/transactions?
    # need to store/buffer?
    # so far, we only need to deal with
    # - end
    # - plot
    # - get input [really?]

    def process(self, action):
        self.post(action)
        # use select to see if there is requests
        return self.retreive()

    def post(self, action):
        self._actions.append(action)
        self.task_queue.put(action)

    def retreive(self):
        result = self.done_queue.get()
        self._actions.pop()
        if isinstance(result, RemoteError):
            if self.debug:
                print(result)
        return result

    def clear(self):
        while len(self._actions) > 0:
            self.retreive()

    def call(self, func, *args, **kwargs):
        action = dict(
            action = 'call',
            func = func,
            args = args,
            kwargs = kwargs,
            )
        return self.process(action)

    def terminate(self, timeout = None):
        # remote here is the kepler process itself, need to adjust for remote
        if not self.remote.is_alive():
            return
        action = dict(
            func = 'terminate',
            )
        self.task_queue.put(action)
        self.task_queue.put('STOP')
        result = self.done_queue.get()
        print(f' [terminate] {result}')
        self.task_queue.join()
        self.remote.join(timeout)
        if self.remote.is_alive():
            self.remote.terminate()
        assert not self.remote.is_alive(), \
               f'Remote process {self.remote.pid} is still running.'
        print(f' [terminate] process {self.remote.name} with PID {self.remote.pid} ended.')

    def execute(self, cmd):
        action = dict(
            action = 'execute',
            cmd = cmd)
        return self.process(action)

    def __getattr__(self, attr):
        get_object = self.remote_arrays ^ attr.endswith('_')
        if attr.endswith('_'):
            attr = attr[:-1]
        if get_object:
            action = dict(
                action = 'type',
                attr = attr,
                )
            rtype = self.process(action)
            if issubclass(rtype, np.ndarray):
                return ProxyArray(self, attr)
        action = dict(
            action = 'get',
            attr = attr,
            )
        return self.process(action)

    def __setattr__(self, attr, value):
        if not attr in self._slots:
            action = dict(
                action = 'set',
                attr = attr,
                value = value,
                )
            return self.process(action)
        return super().__setattr__(attr, value)

    def __delitem__(self, key, value):
        raise NotImplementedError('Deleting remote items not allowed.')

    def __getitem__(self, key):
        action = dict(
            action = 'get',
            )
        if not isinstance(key, tuple):
            attr = key
        else:
            attr = key[0]
            action['index'] = key[1:]
        if not isinstance(attr, str):
            raise TypeError('Require name of remote string.')
        action['attr'] = attr
        val = self.process(action)
        if isinstance(val, RemoteAttributeError):
            raise AttributexError from val
        if isinstance(val, RemoteIndexError):
            raise IndexError from val
        return val

    def __missing__(self, key):
        raise NotImplementedError()

    def __setitem__(self, key, value):
        action = dict(
            action = 'set',
            value = value,
            )
        if not isinstance(key, tuple):
            attr = key
        else:
            attr = key[0]
            action['index'] = key[1:]
        if not isinstance(attr, str):
            raise TypeError('Require name of remote string.')
        action['attr'] = attr
        val = self.process(action)
        if isinstance(val, RemoteAttributeError):
            raise AttributeError from val
        if isinstance(val, RemoteIndexError):
            raise IndexError from val
        return val

    # This should be called cycle
    def s(self, n = 1, wait = True):
        """
        todo - add/update cycle process to do plots
        """
        action = dict(
            action = 'call',
            func = 'cycle',
            kwargs = dict(
                n = n,
                remote = True,
                ),
            )
        if wait:
            return self.process(action)
        self.post(action)

    def __del__(self):
        self.terminate()
        try:
            super().__del__()
        except:
            pass

    @property
    def data(self):
        return ProxyData(self)

    # these are identical with Kepler
    def plot(self, *args, **kwargs):
        self.kp.plot(*args, **kwargs)

    def closewin(self):
        self.kp.closewin()

    def pt(self, mode):
        p = kepplots.get(mode, None)
        if p is not None:
            return p(kd)



class Proxy(ProxyTemplate): # should be LocalProcessProxy
    def setup_process(self, params):
        task_queue = JoinableQueue()
        done_queue = Queue()
        call_queue = JoinableQueue()
        repl_queue = Queue()
        remote = KepProcess(
            task_queue,
            done_queue,
            call_queue,
            repl_queue,
            params,
            )
        remote.start()
        self.remote = remote
        self.task_queue = task_queue
        self.done_queue = done_queue
        self.call_queue = call_queue
        self.repl_queue = repl_queue

class ProxyArray(object):
    """
    interface object to access array on remote object

    Mostly useful on assignment, but also for arrays that remain connected.

    difficult to implement full numpy compatibility

    currently slice operations returns local numpy array only
    """
    _slots = ('proxy', 'name', 'debug')
    def __init__(self, proxy, name, debug = True):
        self.proxy = proxy
        self.name = name
        self.debug = debug
    def __getitem__(self, key):
        if not isinstance(key, tuple):
            key = key,
        arg = (self.name,) + key
        return self.proxy.__getitem__(arg)
    def __setitem__(self, key, value):
        if not isinstance(key, tuple):
            key = key,
        arg = (self.name,) + key
        return self.proxy.__setitem__(arg, value)
    def __getattr__(self, attr):
        if self.debug:
            print(' [' + self.__class__.__name__ + '] ', attr)
        action = dict(
            action = 'get',
            attr = self.name + '.' + attr,
            )
        result =  self.proxy.process(action)
        if isinstance(result, RemoteAttributeError):
            raise AttributeError from result
        return result
    def __setattr__(self, attr, value):
        if attr in self._slots:
            return super().__setattr__(attr, value)
        action = dict(
            action = 'set',
            value = value,
            attr = self.name + '.' + attr,
            )
        return self.proxy.process(action)
    def get(self):
        action = dict(
            action = 'get',
            attr = self.name
            )
        return self.proxy.process(action)
    def __call__(self):
        return self.get()
    def __repr__(self):
        return self.__class__.__name__ + '({})'.format(str(self()))
    def __str__(self):
        return str(self.get())

class ProxyData(DataInterface):
    """
    Proxy object for KeplerData
    """
    _slots = ('proxy',)
    def __init__(self, proxy):
        self.proxy = proxy
    def __getattr__(self, attr):
        action = dict(
            action = 'data',
            attr = attr,
            )
        result =  self.proxy.process(action)
        if isinstance(result, RemoteAttributeError):
            raise AttributeError from result
        return result

def plot(proxy):
    if proxy.__class__.__name__ ==  Proxy.__name__:
        proxy = proxy.data
    return kepplots[1](proxy)
