###############################
# Apply correction factors for the x-sections according to 
# https://indico.cern.ch/event/781231/contributions/3263952/subcontributions/274982/attachments/1795600/2926947/190213_mc_discrep_ppd.pdf
###############################

import yaml
import re

correction_factors = {
	'WJetsToLNu_HT.*(2017|2018)' : {
		'100To200'   : 0.993,
		'200To400'   : 1.002,
		'400To600'   : 1.009,
		'600To800'   : 1.120,
		'800To1200'  : 1.202,
		'1200To2500' : 1.332,
		'2500ToInf'  : 4.200
	},
	'DYJetsToLL.*(2017|2018)' : {
		'100to200'   : 1.000,
		'200to400'   : 0.999,
		'400to600'   : 0.990,
		'600to800'   : 0.975,
		'800to1200'  : 0.907,
		'1200to2500' : 0.833,
		'2500toInf'  : 1.015
	},
	'ZJetsToNuNu.*(2017|2018)' : {
		'100To200'   : 0.994,
		'200To400'   : 0.981,
		'400To600'   : 0.977,
		'600To800'   : 0.975,
		'800To1200'  : 0.916,
		'1200To2500' : 0.880,
		'2500ToInf'  : 1.276
	}	
}

def apply_corrections(xsec_list):
	'''Apply the correction factors to the cross-sections.'''
	datasets = list(xsec_list.keys())
	datasets_corrected = list(correction_factors.keys())
	for dataset in datasets:
		for regex in datasets_corrected:
			if re.match(regex, dataset):
				corrections = correction_factors[regex]
				for ht_bin, correction in corrections.items():
					if ht_bin in dataset:
						print(f'{dataset}: Changing x-section from {xsec_list[dataset]["gen"]} to {xsec_list[dataset]["gen"]*correction:.3f}')
						xsec_list[dataset]['gen'] *= correction

	return xsec_list

def main():
	with open('./xs.yml', 'r') as f:
		xsec_list = yaml.load(f, Loader=yaml.FullLoader)
		updated_xsec_list = apply_corrections(xsec_list)
		
	# Override the original yaml file
	with open('./xs.yml', 'w+') as f:
		yaml.dump(updated_xsec_list, f)

if __name__ == '__main__':
	main()
