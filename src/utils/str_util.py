import codecs
import re


class StrConverter:
    @staticmethod
    def is_chinese(uchar):
        return u'\u4e00' <= uchar <= u'\u9fa5'

    @staticmethod
    def is_number(uchar):
        return u'\u0030' <= uchar <= u'\u0039'

    @staticmethod
    def is_alphabet(uchar):
        return (u'\u0041' <= uchar <= u'\u005a') or (u'\u0061' <= uchar <= u'\u007a')

    @staticmethod
    def file_encode_transform(input_path, input_encode, output_path, output_encode, errors="strict"):
        fw = codecs.open(output_path, 'w', output_encode, errors)
        for line in codecs.open(input_path, 'r', input_encode, errors):
            fw.write(line)
        fw.close()

    @staticmethod
    def extract_file_name(input_name):
        output_name = re.sub(r'\s+', '_', input_name)
        output_name = re.sub(r'[^a-zA-Z0-9_\.]', '', output_name)
        output_name = re.sub(r'\.(txt|.tsv|.csv)', '', output_name)
        if output_name.endswith('result'):
            return output_name[:-6]
        return output_name
