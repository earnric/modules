from cmd import Cmd

class KeplerCmdLoop(Cmd):
    def __init__(self, k, d):
        self.k = k
        self.d = d
        self.prompt = self.d.nameprob + '> '
        super().__init__()

    def default(self, cmdline):
        # filter "s" and "g" commands
        cmd = cmdline.strip().split(' ')
        if cmdline.strip() == 's':
            self.k.s()
            return
        if cmdline.strip() == 'x':
            return True
        if cmd[0] == 's':
            try:
                n = int(cmd[1])
                self.k.cycle(n)
                return
            except:
                pass
        if cmdline.strip() == 'g':
            self.k.g()
            return
        self.k.execute(cmdline)
