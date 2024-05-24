import datetime

def get_time() -> str:
    return datetime.datetime.now().strftime("%d-%b-%Y; %H:%M:%S")
