# -*- coding: utf-8 -*-

class ProgramResults(object):
    """
    This class contains output information from high level programs to easily return a single varible.
    Eventually, logging will be part of this class.
    """
    def __init__(self,name=''):
        self.program_name = name
        self.program_version = None
        self.input_settings = dict() # this is short inputs such as styles, or modes, not whole input vectors. Rule of thumb is the key value pair must fit on one line

    @classmethod
    def new_log(cls):
        return [] # probably make this a class inheriting from list?

    def to_dict(self,names=None):
        if names==None:
            return vars(self)
        else:
            out = {}
            for name in names:
                out[name] = getattr(self, name) # Crash here is operator error
            return out

    def to_list(self,names):
        out = []
        for name in names:
            out.append(getattr(self, name)) # Crash here is operator error
        return out