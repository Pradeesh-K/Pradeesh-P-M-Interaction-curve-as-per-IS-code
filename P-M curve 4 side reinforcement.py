import math
import numpy as np
import pandas as pd
import scipy as sc
import matplotlib
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d
r = [0.05, 0.1, 0.15, 0.2]                                      # the d'/D values
fck = 25                                                        # the grade of concrete used.
st_steel = [0, 0.00144, 0.00163, 0.00192, 0.00241, 0.00276, 0.0038]   # these stress and strain values are for Fe 415
fs_steel = [0, 288.7, 306.7, 324.8, 342.8, 351.8, 360.9]
# st_500 = [0, 0.00174, 0.00195, 0.00226, 0.00277, 0.00312, 0.00417]  # these stress and strain values are for Fe 500
# fs_500 = [0,  347.8, 369.6, 391.3, 413, 423.9, 434.8]
st_pred = interp1d(st_steel, fs_steel, kind = "linear")
c = 40                                                          # the assumed clear cover
p = np.arange(0.02, 0.28, 0.02)                                 # p/fck values ranging from 0.02 to 0.26
d_val = [973.72, 1069.97, 1157.77, 1243.83, 1330.99, 1421.01, 1515.25, 1614.94, 1721.31, 1835.66, 1959.47, 2094.45, 2242.62, 439.18, 457.75, 473.10, 486.86, 499.67, 511.84, 523.57, 534.98, 546.16, 557.18, 568.07, 578.89, 589.66, 283.53, 291.15, 297.29, 302.67, 307.57, 312.14, 316.46, 320.59, 324.57, 328.43, 332.19, 335.86, 339.45, 209.34, 213.47, 216.74, 219.59, 222.16, 224.53, 226.76, 228.87, 230.90, 232.84, 234.72, 236.55, 238.33]
dia_val = [17.37, 27.00, 35.78, 44.38, 53.10, 62.10, 71.52, 81.49, 92.13, 103.57, 115.95, 129.45, 144.26, 7.84, 11.55, 14.62, 17.37, 19.93, 22.37, 24.71, 27.00, 29.23, 31.44, 33.61, 35.78, 37.93, 5.06, 7.34, 9.19, 10.80, 12.27, 13.64, 14.94, 16.18, 17.37, 18.53, 19.66, 20.76, 21.84, 3.74, 5.39, 6.70, 7.84, 8.86, 9.81, 10.70, 11.55, 12.36, 13.14, 13.89, 14.62, 15.33]
bars = 20                                                       # total number of bars assumed
bar_each_side = int(bars/ 4)+1                                  # number of bars in each side
# code to obtain Column Dimension and rebar diameter such that they satisfy both d'/D and p/fck criteria.
# once done the resultant list is extracted and saved. Further code is written to plot data only
# for ratio in r:
#     for percents in p:
#         d = sp.symbols('D')
#         root = sp.solveset(((fck* percents/(bars* 100 * pi)) - ratio**2) * d ** 2 + 2 * ratio * c * d - c ** 2, d )
#         d_val.append(round(max(root.args[0], root.args[1]),2))
#         dia_val.append(round(2 * (ratio * round(max(root.args[0], root.args[1]),2) - c), 2))

c1 = {                                                          # coefficients of area of stress block for Neutral Axis lying outside the column
    1: 0.361,
    1.05 : 0.374,
    1.1 : 0.384,
    1.2: 0.399,
    1.3: 0.409,
    1.4: 0.417,
    1.5: 0.422,
    2.0: 0.435,
    2.5: 0.440,
    3.0: 0.442,
    4.0: 0.444
}
c2 = {                                                           # coefficients of distance of centroid from highly compressed edge for Neutral Axis lying outside the columns
    1.0: 0.416,
    1.05 : 0.432,
    1.1 : 0.443,
    1.2: 0.458,
    1.3: 0.468,
    1.4: 0.475,
    1.5: 0.480,
    2.0: 0.491,
    2.5: 0.495,
    3.0: 0.497,
    4.0: 0.499
}

kv = [0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1,
      1.2, 1.3, 1.4, 1.5, 2.0, 2.5, 3.0, 4.0]        # The ratio of depth of neutral axis to the depth of columnlength of kv  = 29
p_arr = np.ones((52, len(kv)), dtype=float)         # editing this line too p_arr = np.ones((13, len(kv)), dtype=float)
m_arr = np.ones((52, len(kv)), dtype=float)
for i in range(0,52,1):                             # this loop produces normalized P_M Curves for the given values of d'/D
    dia = dia_val[i]
    p_list = []
    m_list = []
    for k in kv:                                    # this loop produces a single normalized P_M Curve
        # print(f' k {k} kd = {k * d_val[i]}')
        if k < 1:
            p_val = 0.36 * k
            m_val = 0.36 * k * (0.5 - 0.416 * k)
            for l in range(bar_each_side):          # iterating through each rebar position to compute the stress and strain for both steel and concrete
                pos = r[i // 13] * d_val[i] + l * (d_val[i] - 2 * r[i// 13] * d_val[i]) / (bar_each_side - 1)
                if k * d_val[i] <= (350 / (350 + 100000 * max(st_steel))) * d_val[i] * (1 - r[i// 13]):
                    # In this case steel reaches yield strain, so concrete doesn't reach 0.0035
                    strain =round(((k * d_val[i] - pos) * max(st_steel) / (d_val[i] - r[i// 13] * d_val[i] - k * d_val[i])), 5)
                    if strain > 0:
                        # print(f' k, d, kd, pos =  {k},   {d_val[i]}, {k * d_val[i]} {pos} strain = {strain}')
                        steel_stress = round(-1 * st_pred(strain), 3)
                    else:
                        # print(f' k, d, kd, pos =  {k},   {d_val[i]}, {k * d_val[i]} {pos} strain = {strain}')
                        steel_stress = round(1 * st_pred(-1 * strain), 3)
                    if pos <= k * d_val[i]:   # if the element lies within N/A
                        if strain >= 0.002: # conc strain > 0.002
                            conc_stress = -0.446 * fck
                        else:
                            conc_stress = round(-0.446 * fck * math.sqrt((k * d_val[i] - pos) * ( 100000 * max(st_steel)/200 )/(d_val[i] * (1 -  r[i // 13] - k))), 3)
                    else:
                        conc_stress = 0
                else:     # now concrete reaches stress
                    # print(f' concrete reaches 0.0035')
                    strain =(k * d_val[i] - pos) * 0.0035 / (k * d_val[i])
                    if strain > 0:
                        steel_stress = round(-1 * st_pred(strain), 3)
                    else:
                        steel_stress = round(1 * st_pred(-1 * strain), 3)
                    if pos <= k * d_val[i]:   # if the element lies within N/A
                        if strain >= 0.002: # conc strain > 0.002
                            conc_stress = -0.446 * fck
                        else:
                            conc_stress = round(-0.5 * fck * math.sqrt(7 * (k * d_val[i] - pos)/(k * d_val[i])), 3)  # made mistake - wrongly put 4/7 kd so got imaginary values
                    else:
                        conc_stress = 0
                if l == 0 or l == 5:
                    p_val = p_val + ((bar_each_side * (3.14 / 4) * dia ** 2) / (fck * d_val[i] ** 2)) * -1 * (steel_stress - conc_stress)           # mistake as -1 is missed.       the 100 is removed in calculations
                    m_val = m_val + (((bar_each_side * (3.14 / 4) * dia ** 2) / (fck * d_val[i] ** 2)) * -1 * (steel_stress - conc_stress)) * ((0.5 * d_val[i] - pos)/ d_val[i])
                else:
                    p_val = p_val + ((2 * (3.14 / 4) * dia ** 2) / (fck * d_val[i] ** 2)) * -1 * (steel_stress - conc_stress)
                    m_val = m_val + (((2 * (3.14 / 4) * dia ** 2) / (fck * d_val[i] ** 2)) * -1 * (steel_stress - conc_stress)) * ((0.5 * d_val[i] - pos) / d_val[i])

        else:
            p_val = c1.get(k)
            m_val = c1.get(k) * (0.5 - c2.get(k))
            for l in range(bar_each_side):
                pos = r[i // 13] * d_val[i] + l * (d_val[i] - 2 * r[i // 13] * d_val[i]) / (bar_each_side - 1)
                strain = (k * d_val[i] - pos) * 0.002 / (k * d_val[i] - (3/7) * d_val[i] )
                steel_stress = round(-1 * st_pred(strain), 3)
                if strain >= 0.002:
                    conc_stress = -0.446 * fck    # i made a mistake here, concrete stress should have been 0.446 for strains > 0.002  .. # mistake - the limiting condition is stress block is 0.0035 at edge with 0 at other edge. However, compressive stress aqt edge decreases with
                else:                               # increase in N/A. Just the fulcrum remains constant
                    conc_stress = round(- 0.446 * fck * math.sqrt((k * d_val[i] - pos)/(k * d_val[i] - (3/7) * d_val[i])), 3)  # mistake: concrete stress, parabolic only upto 0.002 strain , beyond it is straight line. Mistake: missed 4/7 in the multiplier
                if l == 0 or l == 5:
                    p_val = p_val + ((bar_each_side * (3.14 / 4) * dia ** 2) / (fck * d_val[i] ** 2)) * -1 * (steel_stress - conc_stress)
                    m_val = m_val + (((bar_each_side * (3.14 / 4) * dia ** 2) / (fck * d_val[i] ** 2)) * -1 * (steel_stress - conc_stress)) * ((0.5 * d_val[i] - pos) / d_val[i])      # missed -1
                else:
                    p_val = p_val + ((2 * (3.14 / 4) * dia ** 2) / (fck * d_val[i] ** 2)) * -1 * (steel_stress - conc_stress)
                    m_val = m_val + (((2 * (3.14 / 4) * dia ** 2) / (fck * d_val[i] ** 2)) * -1 * (steel_stress - conc_stress)) * ((0.5 * d_val[i] - pos) / d_val[i])             # missed -1
        p_list.append(p_val)
        m_list.append(m_val)
    p_arr[i, :] = p_list           # editing this line and the next, originally p_arr[i % 13, :]
    m_arr[i, :] = m_list
#
pm_arr = np.concatenate((m_arr, p_arr), axis=1)
# Saving the interaction curve data points in excel sheet for column calculations
# df = pd.DataFrame(pm_arr)
# df.to_excel('Fe-500, second part.xlsx')
# print(f'p = {p_arr}')
# print(f'm = {m_arr}')
# plotting the results
fig, ax = plt.subplots(figsize=(12, 8))
plt.ylim(0, 1.2)
for i in range(13):
    ax.plot(m_arr[i, :], p_arr[i, :], marker="x")     # , label=p[i]
plt.title('Compression with bending: Reinforcement on 4 sides, Fy = 415 N/mm2 and d\'/D = 0.2')
# ax.legend()
plt.show()
