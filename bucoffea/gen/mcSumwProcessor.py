import coffea.processor as processor
from awkward import JaggedArray
import numpy as np
from coffea.processor.accumulator import dict_accumulator


def empty_array_100():
    return np.zeros(110)

class mcSumwProcessor(processor.ProcessorABC):
    """
    Processor that stores the nanoaod MC sum of weights for a data set.
    """
    def __init__(self):
        self._acc = dict_accumulator()
        self._acc['sumw'] = processor.defaultdict_accumulator(float)
        self._acc['sumw_scale'] = processor.defaultdict_accumulator(empty_array_100)
        self._acc['sumw_pdf'] = processor.defaultdict_accumulator(empty_array_100)

    @property
    def accumulator(self):
        return self._acc

    def process(self, df):
        out = self.accumulator.identity()
        dataset = df['dataset']
        out['sumw'][dataset] += df['genEventSumw'].sum()

        # pdf_weights = JaggedArray.fromcounts(df['nLHEPdfSumw'], df['LHEPdfSumw'])
        # scale_weights = JaggedArray.fromcounts(df['nLHEScaleSumw'], df['LHEScaleSumw'])

        if 'nLHEPdfSumw' in df:
            out['sumw_pdf'][dataset][:df['nLHEPdfSumw']] += df['LHEPdfSumw']
        if 'nLHEScaleSumw' in df:
            out['sumw_scale'][dataset][:df['nLHEScaleSumw']] += df['LHEScaleSumw']

        return out

    def postprocess(self, acc):
        return acc