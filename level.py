# -*- coding: utf-8 -*-


## /********************************************************/
## /********************************************************/
##                     LEVEL OF THEORY    
## Function: Contains the electronic struture methods (func-
## tionals, MPn, basis set), used in Gaussian 09 and 16.
##
## References:
##
## 1. http://gaussian.com/dft/
## 2. http://gaussian.com/mp/
## 3. http://gaussian.com/basissets/
##
## /********************************************************/
## /********************************************************/


functionals = ["B3LYP", "B3P86", "B3PW91", "O3LYP", "APFD",
               "APF", "LC-wHPBE", "LC-wPBE", "CAM-B3LYP",
               "wB97XD", "wB97", "wB97X", "LC-BLYP", "MN15",
               "M11", "SOGGA11X", "N12SX", "MN12SX", "PW6B95",
               "PW6B95D3", "M08HX", "M06", "M06HF", "M062X",
               "M05", "M052X", "PBE1PBE", "HSEH1PBE", "OHSE2PBE",
               "OHSE1PBE", "PBEh1PBE", "B1B95", "B1LYP", "mPW1PW91",
               "mPW1LYP", "mPW1PBE", "mPW3PBE", "B98", "B971",
               "B972", "TPSSh", "tHCTHhyb", "BMK", "HISSbPBE",
               "X3LYP", "BHandH", "BHandHLYP", "B2PLYP", "B2PLYPD",
               "B2PLYPD3", "mPW2PLYP", "mPW2PLYPD", "DSDPBEP86", "PBE0DH",
               "PBEQIDH", "TPSSTPSS"]

MPn = ["MP2", "MP3", "MP4", "MP4(DQ)", "MP4(SDQ)", "MP4(SDTQ)",
       "MP5", "ROMP2", "ROMP3", "ROMP4"]

CC = ["CC", "CCD", "CCSD", "QCID", "CCSD(T)", "CCSD-T"]

basis_sets = ["STO-3G", "3-21G", "3-21+G", "6-21G", "6-21G*", "6-21G(d,p)",
              "6-21G**", "4-31G", "4-31G*", "4-31G(d,p)", "4-31G**", "6-31G", "6-31+G(d,p)",
              "6-31+G", "6-31++G", "6-31G(d)", "6-31G(d,p)", "6-31G**",
              "6-31G(d')", "6-31G(d',p')", "6-31+G(d'f)", "6-311G", "6-311G(d,p)", "6-311+G", "6-311+G(2d,p)",
              "6-311++G", "6-311++G(d,p)", "6-311++G**", "6-311G(3d,2p)", "D95V", "D95V(d)", "D95V*",
              "D95V(d,p)", "D95V**", "D95", "D95(d)", "D95(d,p)",
              "D95**", "SHC", "SHC*", "CEP-4G", "CEP-4G*", "CEP-31G", "CEP-31G*",
              "CEP-121G", "CEP-121G*", "LanL2MB", "LanL2DZ", "SDD", "SDDAll", "cc-pVDZ",
              "Aug-cc-pVDZ", "cc-pVTZ", "Aug-cc-pVTZ", "cc-pVQZ",  "Aug-cc-pVQZ",
              "cc-pV5Z", "Aug-cc-pV5Z", "cc-pV6Z", "Aug-cc-pV6Z", "Def2SV", "Def2SVP",
              "Def2SVPP", "Def2TZV", "Def2TZVP", "Def2TZVPP", "Def2QZV",
              "Def2QZVP", "Def2QZVPP", "QZVP", "Def2SVPP", "MidiX", "EPR-II",
              "EPR-III", "UGBS", "UGBS1P", "UGBS2V", "MTSmall", "DGDZVP",
              "DGDZVP2", "DGTZVP", "CBSB7"]

