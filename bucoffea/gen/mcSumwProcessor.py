import coffea.processor as processor

class mcSumwProcessor(processor.ProcessorABC):
    """
    Processor that stores the nanoaod MC sum of weights for a data set.
    """
    def __init__(self):
        self._genw = processor.defaultdict_accumulator(float)

    @property
    def accumulator(self):
        return self._genw

    def process(self, df):
        out = self.accumulator.identity()
        out[df['dataset']] += df['genEventSumw'].sum()
        return out

    def postprocess(self, acc):
        return acc