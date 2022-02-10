
import numpy as np


def cal_to_julian(yr, mon, day, hr, min, sec):

	JD = 367 * yr - int((7 * (yr + int((mon + 9) / 12))) / 4) + int(275 * mon / 9) + day + 1721013.5 + ((hr + (min + (sec / 60)) / 60) / 24)

	return JD


def MJD(JD):

	return JD - 2400000.5


def julian_to_cal(JD):

	#Algorithm valid 1900 - 2100

	LMonth = np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])

	T_1900 = (JD - 2415019.5) / 365.25
	year = 1900 + np.trunc(T_1900)

	leap_years = np.trunc((year - 1900 - 1) * 0.25)
	if year % 4 == 0:
		LMonth[2 - 1] = 29

	days = (JD - 2415019.5) - ((year - 1900) * 365.0 + leap_years)
	if days < 1.0:
		year = year - 1
		leap_years = np.trunc((year - 1900 - 1) * 0.25)
		days = (JD - 2415019.5) - ((year - 1900) * 365.0 + leap_years)

	DOY = np.trunc(days)

	mon = 0
	day = 0
	sum_days = 0
	for i in range(len(LMonth)):

		sum_days += LMonth[i]

		if sum_days > DOY:

			day = DOY - (sum_days - LMonth[i])
			mon = i + 1

			break

	t = (days - DOY) * 24
	hrs = np.trunc(t)
	mins = np.trunc((t - hrs) * 60)
	sec = (t - hrs - (mins / 60)) * 3600

	return [int(year), int(mon), int(day), int(hrs), int(mins), sec]


