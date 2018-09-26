#! /usr/bin/python3.6
import math

def decdeg2hexdeg(decdeg):
    ''' Decimal degrees to hexicondal degrees.
        
        Args:
            decdeg (float): angle in decimal degrees
        Returns:
            tuple (int, int, float): The elements are:
                # 0 -> integer degrees
                # 1 -> integer minutes
                # 2 -> float seconds
    '''
    dd = decdeg
    dd1  = abs(float(dd))
    cdeg = int(dd1)
    minsec = dd1 - cdeg
    cmin = int(minsec * 60)
    csec = (minsec % 60) / float(3600)
    if dd < 0e0: cdeg = cdeg * -1
    return cdeg,cmin,csec

def rad2hexdeg(rad):
    decdeg = math.degrees(rad)
    return decdeg2hexdeg(decdeg)

def normalize(num, lower=0e0, upper=360e0):
    ''' Normalize number in the range [lower, upper)
    '''
    if lower >= upper:
        raise ValueError("Invalid lower and upper limits: (%s, %s)")
    res = num
    if num > upper or num == lower:
        num = lower + abs(num + upper) % (abs(lower) + abs(upper))
    if num < lower or num == upper:
        num = upper - abs(num - lower) % (abs(lower) + abs(upper))

    res = lower if res == upper else num
    return res

if __name__ == "__main__":
    from angles import normalize as nrm
    angle =-3e0
    while angle < math.pi*4.5e0:
        a2 = normalize(angle, 0e0, 2e0*math.pi)
        a1 = nrm(angle, 0e0, 2e0*math.pi)
        if abs(a1-a2) > 1e-10:
            print("Different result for angle {:+6.2f}; normalize returns {:+6.2f}, normalize2 returns {:+6.2f}".format(math.degrees(angle), math.degrees(a1), math.degrees(a2)))
        angle += .003
    angle =-3e0
    while angle < math.pi*4.5e0:
        a2 = normalize(angle, -math.pi, math.pi)
        a1 = nrm(angle, -math.pi, math.pi)
        if abs(a1-a2) > 1e-10:
            print("Different result for angle {:+6.2f}; normalize returns {:+6.2f}, normalize2 returns {:+6.2f}".format(math.degrees(angle), math.degrees(a1), math.degrees(a2)))
        angle += .003
    print("All done")
