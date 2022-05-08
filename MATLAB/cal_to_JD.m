
function JD = cal_to_JD(year, month, day, hrs, mins, sec)

    JD = 367 * year - fix((7 * (year + fix((month + 9) / 12))) / 4) + fix(275 * month / 9) + day + 1721013.5 + ((hrs + (mins + (sec / 60)) / 60) / 24);

end