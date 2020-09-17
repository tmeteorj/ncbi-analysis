import datetime


def get_time():
    return datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
