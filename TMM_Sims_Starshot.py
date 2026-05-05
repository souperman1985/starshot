import numpy as np
import pandas as pd
from scipy.optimize import least_squares
from scipy.interpolate import interp1d
import math

def Stack(lamb, UnitLayer, thicknesses, Substrate, SubstrateThickness, NumLayers=1, phi=0, n=None, TM=0,
          output='Absorptance', n_inc=1, SubstrateAbsorption=True):
    t = np.hstack([thicknesses]*NumLayers)
    t = np.hstack((np.array([0]), t, SubstrateThickness))

    lamb = lamb

    if not isinstance(n, np.ndarray):
        #Load in index data
        specs = pd.read_excel('MatData.xlsx', skiprows=2, header=None).to_numpy()
        names = pd.read_excel('MatData.xlsx', nrows=2, header=None).to_numpy()

        #Construct Unit Layer
        n = np.zeros((len(UnitLayer), np.shape(lamb)[0])).astype(np.complex128)
        for i, mat in enumerate(UnitLayer):
            c = np.argwhere(names[0, :] == mat)[0][0]
            f1 = interp1d(specs[:, c-1], specs[:, c])
            f2 = interp1d(specs[:, c-1], specs[:, c+1])
            n[i, :] = f1(lamb) + 1j*f2(lamb)

        #Stack unit layers
        n = np.vstack([n]*NumLayers)

        #Sandwich between air and substrate
        if not isinstance(n_inc, np.ndarray):
            n = np.vstack((np.ones((1, np.shape(lamb)[0]))*n_inc, n)) #put air on the top
        else:
            n = np.vstack((n_inc, n))

        for sub in Substrate:
            mat = np.argwhere(names[0, :] == sub)[0][0]
            f1 = interp1d(specs[:, mat-1], specs[:, mat])
            f2 = interp1d(specs[:, mat-1], specs[:, mat+1])
            n = np.vstack((n, f1(lamb) + 1j*f2(lamb)))

    q = (n**2 - np.sin(phi))**0.5

    #Calculate transfer matrices, and field at each wavelength and position
    T = []
    R = []
    if TM == 0:
        for l, wl in enumerate(lamb):
            # Calculate transfer matrices for incoherent reflection and transmission at the first interface
            r = (q[0, l] - q[1, l])/(q[0, l] + q[1, l]) #Eq. 2a
            m = (2*q[0, l])/(q[0, l] + q[1, l]) #Eq. 2b
            S = np.array([[1, r], [r, 1]])/m #First term of Eq. 8

            #Go thru the terms in Eq. 8 and update S
            for matindex in range(1, np.shape(t)[0]-1):
                xi = 2*np.pi*q[matindex, l]/wl #Eq. 6
                L_mat = np.array([[np.exp(-1j*xi*t[matindex].astype(np.complex128)), 0], [0, np.exp(1j*xi*t[matindex].astype(np.complex128))]]) #Eq. 5
                r = (q[matindex, l] - q[matindex + 1, l])/(q[matindex, l] + q[matindex + 1, l]) #Eq. 2a
                m = (2*q[matindex, l])/(q[matindex, l] + q[matindex + 1, l]) #Eq. 2b
                I_mat = np.array([[1, r], [r, 1]])/m #Eq. 1
                S = np.matmul(S, L_mat) #Update S
                S = np.matmul(S, I_mat)

            #Added by Jason
            matindex += 1
            xi = 2*np.pi*q[matindex, l]/wl
            L_mat = np.array([[np.exp(-1j*xi*t[matindex].astype(np.complex128)), 0], [0, np.exp(1j*xi*t[matindex].astype(np.complex128))]]) #Eq. 5
            S = np.matmul(S, L_mat)
            #End of added by Jason

            R.append(S[1, 0]/S[0, 0])
            T.append(1/S[0, 0])
    else:
        for l, wl in enumerate(lamb):
            # Calculate transfer matrices for incoherent reflection and transmission at the first interface
            r = ((q[0, l]*(n[1, l]**2)) - ((n[0, l]**2)*q[1, l]))/((q[0, l]*(n[1, l]**2)) + ((n[0, l]**2)*q[1, l])) #Eq. 2a
            m = (2*q[0, l]*n[0, l]*n[1, l])/((q[0, l]*(n[1, l]**2)) + ((n[0, l]**2)*q[1, l])) #Eq. 2b
            S = np.array([[1, r], [r, 1]])/m #First term of Eq. 8

            #Go thru the terms in Eq. 8 and update S
            for matindex in range(1, np.shape(t)[0]-1):
                xi = 2*np.pi*q[matindex, l]/wl #Eq. 6
                L_mat = np.array([[np.exp(-1j*xi*t[matindex]), 0], [0, np.exp(1j*xi*t[matindex])]]) #Eq. 5
                r = ((q[matindex, l]*(n[matindex+1, l]**2)) - ((n[matindex, l]**2)*q[matindex + 1, l]))/((q[matindex, l]*(n[matindex+1, l]**2)) + ((n[matindex, l]**2)*q[matindex + 1, l])) #Eq. 2a
                m = (2*q[matindex, l]*n[matindex, l]*n[matindex+1, l])/((q[matindex, l]*(n[matindex+1, l]**2)) + ((n[matindex, l]**2)*q[matindex + 1, l])) #Eq. 2b
                I_mat = np.array([[1, r], [r, 1]])/m #Eq. 1
                S = np.matmul(S, L_mat) #Update S
                S = np.matmul(S, I_mat)

            if SubstrateAbsorption == True:
                matindex += 1
                xi = 2*np.pi*q[matindex, l]/wl
                L_mat = np.array([[np.exp(-1j*xi*t[matindex].astype(np.complex128)), 0], [0, np.exp(1j*xi*t[matindex].astype(np.complex128))]]) #Eq. 5
                S = np.matmul(S, L_mat)

            R.append(S[1, 0]/S[0, 0])
            T.append(1/S[0, 0])

    R = np.array(R)
    T = np.array(T)
    Abs = 1 - (np.abs(R)**2)*n_inc - (np.abs(T)**2)*np.real(n[-1, :])

    if output == 'Absorptance':
        Result = Abs
    elif output == 'Reflectance':
        Result = np.abs(R)**2
    elif output == 'Reflection':
        Result = R
    elif output == 'Transmittance':
        Result = (np.abs(T)**2)*np.real(n[-1, :])
    elif output == 'Transmission':
        Result = T
    elif output == 'R and A':
        Result = np.vstack((np.abs(R)**2, Abs))
    elif output == 'All':
        Result = np.vstack((Abs, np.abs(R)**2, np.abs(T)**2)*np.real(n[-1, :]))
    else:
        print('Invalid output input')

    return Result