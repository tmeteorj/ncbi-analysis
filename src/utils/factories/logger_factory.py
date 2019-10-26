import datetime
import time


class LoggerFactory:

    def __init__(self, time_distance=5):
        self.time_distance = time_distance
        self.start_time = 0
        self.last_time = None

    def info_with_expire_time(self, msg, solve, total):
        if solve == 0 or self.start_time == 0:
            self.last_time = None
            self.start_time = time.time()
        if not self.last_time or time.time() - self.last_time >= self.time_distance:
            remain = None
            if solve > 0 and total > 1:
                remain = (time.time() - self.start_time) / solve * (total - solve)
                if remain < 0:
                    remain = 0
            if remain is None:
                print(self.get_time() + ": " + msg)
            else:
                print("%s: %s, remain %02d:%02d:%02d" % (
                    self.get_time(), msg, int(remain / 3600), int(remain / 60 % 60), int(remain % 60)))
            self.last_time = time.time()

    def info_per_time(self, msg):
        if self.start_time == 0:
            self.last_time = None
            self.start_time = time.time()
        if not self.last_time or time.time() - self.last_time >= self.time_distance:
            print(self.get_time() + ": " + msg)
            self.last_time = time.time()

    @staticmethod
    def info(msg):
        print(LoggerFactory.get_time() + ": " + msg)

    @staticmethod
    def get_time():
        return datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
