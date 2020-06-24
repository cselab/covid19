import datetime

def date_fromisoformat(s):
    """Parse a YYYY-MM-DD date format to a datetime.date object.

    Reproduces Python 3.7's datetime.date.fromisoformat.
    """

    dstr = s[0:10] # sometimes s comes with hr, sec and ms
    year, month, day = dstr.split('-')
    return datetime.date(int(year), int(month), int(day))
