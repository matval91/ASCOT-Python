import numpy as np

def cocos_transform(eqd, cocosin_f, cocosout_f, sigma_ip_out=1, sigma_b0_out=1):
    """ cocos transformations
    This function converts the magnetic input from their starting cocos to cocos 3 (needed by ascot5)

    Parameters:
        COCOS (int): input cocos. Now useable only 2,3,4,7,12,13,14,17
    Attributes:
        None

    """

    print("COCOS tranformation from "+str(cocosin_f)+" to "+str(cocosout_f))
    cocos_keys = ['sigma_Bp', 'sigma_RphiZ', 'sigma_rhothetaphi',\
                  'sign_q_pos', 'sign_pprime_pos', 'exp_Bp']
    pi = np.pi
    
    eqdout = eqd

    cocosin  = dict.fromkeys(cocos_keys)
    cocosout = dict.fromkeys(cocos_keys)
    cocosin['exp_Bp'] = 0 # if cocosin>=10, this should be 1
    if cocosin_f>10:
        cocosin['exp_Bp'] = +1 # if cocosin>=10, this should be 1


    if cocosin_f==2 or cocosin_f==12:
        #These cocos are for CHEASE (2) - IN
        cocosin['sigma_Bp'] = +1
        cocosin['sigma_RphiZ'] = -1
        cocosin['sigma_rhothetaphi'] = +1
        cocosin['sign_q_pos'] = +1
        cocosin['sign_pprime_pos'] = -1
    elif cocosin_f==3 or cocosin_f==13:
        #These cocos are for EFIT (3) - IN
        cocosin['sigma_Bp'] = -1
        cocosin['sigma_RphiZ'] = +1
        cocosin['sigma_rhothetaphi'] = -1
        cocosin['sign_q_pos'] = -1
        cocosin['sign_pprime_pos'] = +1
    elif cocosin_f==4 or cocosin_f==14:
        cocosin['sigma_Bp'] = -1
        cocosin['sigma_RphiZ'] = -1
        cocosin['sigma_rhothetaphi'] = -1
        cocosin['sign_q_pos'] = -1
        cocosin['sign_pprime_pos'] = +1
    elif cocosin_f==5 or cocosin_f==15:
        cocosin['sigma_Bp'] = 1
        cocosin['sigma_RphiZ'] = +1
        cocosin['sigma_rhothetaphi'] = +1
        cocosin['sign_q_pos'] = +1
        cocosin['sign_pprime_pos'] = -1  
    elif cocosin_f==7 or cocosin_f==17:
        #These cocos are for LIUQE(17) - IN
        cocosin['sigma_Bp'] = -1
        cocosin['sigma_RphiZ'] = +1
        cocosin['sigma_rhothetaphi'] = +1
        cocosin['sign_q_pos'] = +1
        cocosin['sign_pprime_pos'] = +1
    else:
        print(str(cocosin_f)+" Not Implemented \n")
        raise ValueError

    cocosin['sigma_ip'] = np.sign(eqd.Ip)
    cocosin['sigma_b0'] = np.sign(eqd.B0EXP)

    cocosout['exp_Bp'] = 0 # if cocosin>=10, this should be 1
    if cocosout_f>10:
        cocosout['exp_Bp'] = +1 # if cocosin>=10, this should be 1

    #Defining output cocos
    if cocosout_f==2 or cocosout_f==12:
        #These cocos are for CHEASE (2) - IN
        cocosout['sigma_Bp'] = +1
        cocosout['sigma_RphiZ'] = -1
        cocosout['sigma_rhothetaphi'] = +1
        cocosout['sign_q_pos'] = +1
        cocosout['sign_pprime_pos'] = -1
    elif cocosout_f==3 or cocosout_f==13:
        #These cocos are for EFIT (3) - IN
        cocosout['sigma_Bp'] = -1
        cocosout['sigma_RphiZ'] = +1
        cocosout['sigma_rhothetaphi'] = -1
        cocosout['sign_q_pos'] = -1
        cocosout['sign_pprime_pos'] = +1
    elif cocosout_f==4 or cocosout_f==14:
        cocosout['sigma_Bp'] = -1
        cocosout['sigma_RphiZ'] = -1
        cocosout['sigma_rhothetaphi'] = -1
        cocosout['sign_q_pos'] = -1
        cocosout['sign_pprime_pos'] = +1
    elif cocosout_f==5 or cocosout_f==15:
        cocosout['sigma_Bp'] = 1
        cocosout['sigma_RphiZ'] = +1
        cocosout['sigma_rhothetaphi'] = +1
        cocosout['sign_q_pos'] = +1
        cocosout['sign_pprime_pos'] = -1      
    elif cocosout_f==7 or cocosout_f==17:
        #These cocos are for LIUQE(17) - IN
        cocosout['sigma_Bp'] = -1
        cocosout['sigma_RphiZ'] = +1
        cocosout['sigma_rhothetaphi'] = +1
        cocosout['sign_q_pos'] = +1
        cocosout['sign_pprime_pos'] = +1
    else:
        print(str(cocosout_f)+" Not Implemented \n")
        raise ValueError
        
    cocosout['sigma_ip'] = sigma_ip_out
    cocosout['sigma_b0'] = sigma_b0_out

    # Define effective variables: sigma_Ip_eff, sigma_B0_eff, sigma_Bp_eff, exp_Bp_eff as in Appendix C
    #sigma_Ip_eff = cocosin['sigma_RphiZ'] * cocosout['sigma_RphiZ']
    #sigma_B0_eff = cocosin['sigma_RphiZ'] * cocosout['sigma_RphiZ']
    # Since we want sigmaip and sigmab0 defined, we must use
    sigma_Ip_eff = cocosin['sigma_ip']*cocosout['sigma_ip']
    sigma_B0_eff = cocosin['sigma_b0']*cocosout['sigma_b0']
    sigma_Bp_eff = cocosin['sigma_Bp'] * cocosout['sigma_Bp']
    exp_Bp_eff = cocosout['exp_Bp'] - cocosin['exp_Bp']
    sigma_rhothetaphi_eff = cocosin['sigma_rhothetaphi'] * cocosout['sigma_rhothetaphi']
    # Define input
    F_in = eqd.T
    FFprime_in = eqd.TTprime
    pprime_in = eqd.pprime
    psirz_in = eqd.psi
    psiaxis_in = eqd.psiaxis
    psiedge_in = eqd.psiedge
    q_in = eqd.q
    b0_in = eqd.B0EXP
    ip_in = eqd.Ip
    psigrid_in = eqd.psi_grid
    # Transform
    F = F_in * sigma_B0_eff
    FFprime = FFprime_in*sigma_Ip_eff*sigma_Bp_eff/(2*pi)**exp_Bp_eff
    pprime = pprime_in * sigma_Ip_eff * sigma_Bp_eff / (2*pi)**exp_Bp_eff
    _fact_psi = sigma_Ip_eff * sigma_Bp_eff * (2*pi)**exp_Bp_eff
    psirz = psirz_in * _fact_psi 
    psigrid = psigrid_in * _fact_psi
    psiaxis = psiaxis_in * _fact_psi
    psiedge = psiedge_in * _fact_psi
    q = q_in * sigma_Ip_eff * sigma_B0_eff * sigma_rhothetaphi_eff
    b0 = b0_in * sigma_B0_eff
    ip = ip_in * sigma_Ip_eff
    # Define output
    eqdout.T = F
    eqdout.TTprime = FFprime
    eqdout.pprime = pprime
    eqdout.psi = psirz
    eqdout.psigrid = psigrid
    eqdout.psiaxis = psiaxis
    eqdout.psiedge = psiedge
    eqdout.q = q
    eqdout.B0EXP = b0
    eqdout.Ip = ip
    return eqdout
