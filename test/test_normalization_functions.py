import unittest
import LongPass.normalization_functions as normalization_functions

class TestNormalizationFunction(unittest.TestCase):

    def load_lrtsp_count(self,filepath):
        counts = []
        with open(filepath,'r') as lrtspfile:
            for line in lrtspfile:
                chrname = line.split('\t')[0]
                pos = line.split('\t')[1]
                count = line.split('\t')[3]
                counts.append([chrname,int(pos),int(count)])

        return sorted(counts, key=lambda x: x[0])

    def load_testnormalizedtpms(self,filepath):
        testnormalizedtpms = []
        with open(filepath,'r') as testnormalizedtpmsfile:
            for testnormalized in testnormalizedtpmsfile:
                testchrname = testnormalized.split('\t')[0]
                testpos = testnormalized.split('\t')[1]
                testnormalizedtpm = round(float(testnormalized.split('\t')[3]),5)
                testnormalizedtpms.append([testchrname,int(testpos),testnormalizedtpm])

        return sorted(testnormalizedtpms,key=lambda x: x[0])   #why need sort? cause the loading of the lrtsp file will change the order of chr name in CAGEr PACKAGE.

    def test_simpletpm(self):
        counts = self.load_lrtsp_count("/home/hongyanhong/TCfinder/test/start0909.lrtsp")
        total_count = sum([item[2] for item in counts])
        normalized_counts = [normalization_functions.simpleTpm(item[2], total_count) for item in counts]
        normalized_list = []
        for count,normalized_count in zip(counts,normalized_counts):
            normalized_list.append([count[0],count[1],normalized_count])
        testnormalizedtpms = self.load_testnormalizedtpms("/home/hongyanhong/TCfinder/test/normalizedTpm.txt")

        self.assertEqual(normalized_list, testnormalizedtpms)

    def test_fit_power_law_to_values(self):
        counts = self.load_lrtsp_count("/home/hongyanhong/TCfinder/test/start0909.lrtsp")
        counts = [count[2] for count in counts]
        slope,intercept = normalization_functions.fit_power_law_to_values(counts)
        with open('/home/hongyanhong/TCfinder/test/test_readscounts_results.txt') as test_readscounts_results_file:
            for line in test_readscounts_results_file:
                sloperesult = float((line.strip('\n').split(',')[0]))
                interceptresult = float(line.strip('\n').split(',')[1])
        self.assertEqual([sloperesult,interceptresult], [round(slope,6),round(intercept,6)])

    def test_fit_values_to_referent_power_law(self):
        chr_pos_counts = self.load_lrtsp_count("/home/hongyanhong/TCfinder/test/start0909.lrtsp")
        counts = [count[2] for count in chr_pos_counts]
        slope,intercept = normalization_functions.fit_power_law_to_values(counts)
        values_fitted = normalization_functions.fit_values_to_referent_power_law(counts,[slope,intercept])
        values_test_results = [[chr_pos_count[0],chr_pos_count[1],value_fitted] for chr_pos_count,value_fitted in zip(chr_pos_counts,values_fitted)]
        values_test_results = sorted(values_test_results, key=lambda x:x[0], reverse=False)
        values_real_results = []
        with open('/home/hongyanhong/TCfinder/test/test_fit_values_to_referent_power_law_results.txt') as file:
            for line in file:
                chro = line.strip('\n').split('\t')[0]
                pos = int(line.strip('\n').split('\t')[1])
                fittedvalue = round(float(line.strip('\n').split('\t')[2]),5)
                values_real_results.append([chro,pos,fittedvalue])
        
        values_real_results = sorted(values_real_results, key=lambda x:x[0], reverse=False)

        self.assertEqual(values_test_results, values_real_results)