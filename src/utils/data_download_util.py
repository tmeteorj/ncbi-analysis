import traceback
from urllib import request

from experiment_config import *


class DataDownloadTool:

    @staticmethod
    def format_data(line):
        key, inter = line[1:].strip().split('/')
        inter = list(map(int, inter.split('-')))
        return key, inter

    @staticmethod
    def get_response(key):
        url = URL_LIB_PREFIX + key
        return request.urlopen(url, timeout=60)

    @staticmethod
    def read_data(uid, download_size, timeout):
        url = URL_DATA_DOWNLOAD % (uid, download_size)
        return request.urlopen(url, timeout=timeout)

    @staticmethod
    def read_uid(res):
        info = res.readlines()[6].decode('utf8').strip().split('/>')[:-1]
        dic = {}
        for data in info:
            meta, name, content = data.strip().split(' ')
            name, content = [x.strip().split('=')[1].strip('\"') for x in [name, content]]
            dic[name] = content
        return dic[KEY_LIB_UID]

    @staticmethod
    def download_data(key, filename):
        uid = None
        for try_time in range(1, MAX_DOWNLOAD_RETRY_TIME + 1):
            try:
                uid = DataDownloadTool.read_uid(DataDownloadTool.get_response(key))
                if not uid:
                    continue
                break
            except:
                print("get uid failed, retry %d/%d, msg : " % (try_time, MAX_DOWNLOAD_RETRY_TIME))
                traceback.print_exc()
        download_size = 1000000 * 100
        timeout = 60
        flag = False
        if uid is None:
            print('get uid failed finally')
            return flag
        for try_time in range(1, MAX_DOWNLOAD_RETRY_TIME + 1):
            try:
                data = DataDownloadTool.read_data(uid, download_size, timeout)
                line_idx = 0
                with open(filename, 'wb') as f:
                    for x in data.readlines():
                        f.write(x)
                        line_idx += 1
                if line_idx < 10:
                    print('download key %s failed for download_size=%dM, timeout=%ds, retry %d/%d' % (
                        key, download_size / 1000000, timeout, try_time, MAX_DOWNLOAD_RETRY_TIME))
                    download_size += 1000000 * 100
                    timeout += 300
                else:
                    flag = True
                    break
            except:
                traceback.print_exc()
        return flag

