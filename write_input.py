# -*- coding: utf-8 -*-

## /********************************************************/
## /********************************************************/
##                     WRITE INPUT 
## Function: Subroutine that writes input files for computa-
## tions (single point, optimization, transition state opti-
## mization, frequency analysis, intrinsic reaction coordi-
## nate (IRC), and volume) with Gaussian (09 and 16).
## /********************************************************/
## /********************************************************/

#-----------------------------------------------------------
#Python modules
import os
from os import path
#Eyring modules
from gaussian_parser import get_xyz_gaussian, get_charge, get_multiplicity
from constants import *
from sort_by_rxcoord import *
from utils import get_output, get_xyz_file
from xyzparser import read_xyz
from get_xyz_gaussian import *
from get_xyz_gaussian_scan import *
#-----------------------------------------------------------
def write_input(mem, keys, charge, mult, atom, x, y, z,
                dirname, outfile, unrestricted = "false"):

    if not os.path.isdir("./" + dirname):
        os.mkdir(dirname)
        
    for i in range(len(atom)):
        basename = os.path.splitext(os.path.basename(outfile[i]))[0]
        if "irc" in keys.lower():
            if "forward" in keys.lower():
                inpname = str(basename) + "_ircfor" + ".inp"
            elif "reverse" in keys.lower():
                inpname = str(basename) + "_ircrev" + ".inp"
            else:
                inpname = str(basename) + "_irc" + ".inp"       
        else:
            inpname = str(basename) + ".inp"
        inpfile = os.path.join(dirname, inpname)
        with open(inpfile,'w') as inp:
            inp.write("%%MEM=%-6s \n" % mem)
            if unrestricted == "false":
                inp.write("#%-6s \n" % keys)
            if unrestricted == "true":
                if int(mult[i]) == 1:
                    inp.write("#%-6s \n" % keys)
                if int(mult[i]) > 1:
                    inp.write("#U%-6s \n" % keys)
            inp.write("\n")
            inp.write("comment\n")
            inp.write("\n")
            inp.write("%-2.0f %2.0f \n" % (charge[i], mult[i]))
            for k in range(len(atom[i][-1])):
                inp.write("%-4s % -6.9f   % -6.9f   % -6.9f     \n"
                          % (atom[i][-1][k], x[i][-1][k], y[i][-1][k], z[i][-1][k]))
            inp.write("\n")
#-----------------------------------------------------------
def write_input_for_xyz(mem, keys, charge, mult, atom, x, y, z,
                dirname, outfile, unrestricted = "false"):

    if not os.path.isdir("./" + dirname):
        os.mkdir(dirname)
        
    for i in range(len(atom)):
        basename = os.path.splitext(os.path.basename(outfile[i]))[0]
        if "irc" in keys.lower():
            if "forward" in keys.lower():
                inpname = str(basename) + "_ircfor" + ".inp"
            elif "reverse" in keys.lower():
                inpname = str(basename) + "_ircrev" + ".inp"
            else:
                inpname = str(basename) + "_irc" + ".inp"       
        else:
            inpname = str(basename) + ".inp"
        inpfile = os.path.join(dirname, inpname)
        with open(inpfile,'w') as inp:
            inp.write("%%MEM=%-6s \n" % mem)
            if unrestricted == "false":
                inp.write("#%-6s \n" % keys)
            if unrestricted == "true":
                if int(mult[i]) == 1:
                    inp.write("#%-6s \n" % keys)
                if int(mult[i]) > 1:
                    inp.write("#U%-6s \n" % keys)
            inp.write("\n")
            inp.write("comment\n")
            inp.write("\n")
            inp.write("%-2.0f %2.0f \n" % (charge, mult))
            for k in range(len(atom[i])):
                inp.write("%-4s % -6.9f   % -6.9f   % -6.9f     \n"
                          % (atom[i][k], x[i][k], y[i][k], z[i][k]))
            inp.write("\n")
#-----------------------------------------------------------            
def write_input_from_irc(mem, keys, charge, mult, atom, x, y, z,
                         dirname, rxcoord, energy, unrestricted = "false",
                         position = None):

    if position is None:
        position = [ ]

    if not os.path.isdir("./" + dirname):
        os.mkdir(dirname)

    for i in range(len(atom)):
        if len(position) == 0:
            num = str(i + 1).zfill(4)
        else:
            num = str(position[i]).zfill(4)
        inpname = str(num) + ".inp"
        inpfile = os.path.join(dirname, inpname)
        with open(inpfile,'w') as inp:
            inp.write("%%MEM=%-6s \n" % mem)
            if unrestricted == "false":
                inp.write("#%-6s \n" % keys)
            if unrestricted == "true":
                if int(mult) == 1:
                    inp.write("#%-6s \n" % keys)
                if int(mult) > 1:
                    inp.write("#U%-6s \n" % keys)
            inp.write("\n")
            inp.write("RxCoord=%-6.5f   DE=%-6.2f \n"
                      % (rxcoord[i], ((energy[i] * HARTREES_TO_CALORIES_MOL) / 1000.0)))
            inp.write("\n")
            inp.write("%-2.0f %2.0f \n" % (charge, mult))
            for j in range(len(atom[i])):
                inp.write("%-4s % -6.9f   % -6.9f   % -6.9f \n"
                          % (atom[i][j], x[i][j], y[i][j], z[i][j]))
            inp.write("\n")
#-----------------------------------------------------------
def write_inputchk_from_irc(mem, keys, charge, mult, atom, x, y, z,
                         dirname, rxcoord, energy, unrestricted = "false",
                         position = None):

    if position is None:
        position = [ ]

    if not os.path.isdir("./" + dirname):
        os.mkdir(dirname)

    CHKdirname = dirname + str("/CHKfiles")

    os.mkdir(CHKdirname)

    for i in range(len(atom)):
        if len(position) == 0:
            num = str(i + 1).zfill(4)
        else:
            num = str(position[i]).zfill(4)

        chkfile = "CHKfiles/" + num + ".chk"

        inpname = str(num) + ".inp"
        inpfile = os.path.join(dirname, inpname)
        with open(inpfile,'w') as inp:
            inp.write("%%MEM=%-6s \n" % mem)
            inp.write("%%chk=%-6s \n" % chkfile)
            if unrestricted == "false":
                inp.write("#%-6s \n" % keys)
            if unrestricted == "true":
                if int(mult) == 1:
                    inp.write("#%-6s \n" % keys)
                if int(mult) > 1:
                    inp.write("#U%-6s \n" % keys)
            inp.write("\n")
            inp.write("RxCoord=%-6.5f   DE=%-6.2f \n"
                      % (rxcoord[i], ((energy[i] * HARTREES_TO_CALORIES_MOL) / 1000.0)))
            inp.write("\n")
            inp.write("%-2.0f %2.0f \n" % (charge, mult))
            for j in range(len(atom[i])):
                inp.write("%-4s % -6.9f   % -6.9f   % -6.9f \n"
                          % (atom[i][j], x[i][j], y[i][j], z[i][j]))
            inp.write("\n")
#-----------------------------------------------------------
def write_input_from_scan(mem, keys, charge, mult, atom, x, y, z,
                         dirname, unrestricted = "false", position = None):

    if position is None:
        position = [ ]

    if not os.path.isdir("./" + dirname):
        os.mkdir(dirname)

    for i in range(len(atom)):
        if len(position) == 0:
            num = str(i + 1).zfill(4)
        else:
            num = str(position[i]).zfill(4)
        inpname = str(num) + ".inp"
        inpfile = os.path.join(dirname, inpname)
        with open(inpfile,'w') as inp:
            inp.write("%%MEM=%-6s \n" % mem)
            if unrestricted == "false":
                inp.write("#%-6s \n" % keys)
            if unrestricted == "true":
                if int(mult) == 1:
                    inp.write("#%-6s \n" % keys)
                if int(mult) > 1:
                    inp.write("#U%-6s \n" % keys)
            inp.write("\n")
            inp.write("comment \n")
            inp.write("\n")
            inp.write("%-2.0f %2.0f \n" % (charge, mult))
            for j in range(len(atom[i])):
                inp.write("%-4s % -6.9f   % -6.9f   % -6.9f \n"
                          % (atom[i][j], x[i][j], y[i][j], z[i][j]))
            inp.write("\n")
#-----------------------------------------------------------
def write_inputs():
    args = sys.argv[1:]
    modulename = os.path.splitext(os.path.basename(sys.argv[0]))[0]
    cwd = os.path.basename(os.getcwd())

    #Parsing script options
    if len(args) >= 6:
        arg1 = str(sys.argv[2])
        dirname = str(sys.argv[3])
        mem = str(sys.argv[4])
        keys = str(sys.argv[5])
        unrestricted = str(sys.argv[6]).lower()
        if arg1 == "g":
            opt = "1"
        elif arg1 == "i":
            opt = "2"
        elif arg1 == "i1":
            opt = "3"
        elif arg1 == "i2":
            opt = "4"
        elif arg1 == "x":
            opt = "5"
        elif arg1 == "s":
            opt = "6"
        elif arg1 == "s1":
            opt = "7"
        elif arg1 == "s2":
            opt = "8"
        elif arg1 == "c":
            opt = "9"
        else:
            print("Error: Invalid option!")
            sys.exit(1)
    else:
        print("Error: Invalid option!") 
        sys.exit(1)

    if opt == '1':
        if len(args) == 6:
            out = sorted(get_output(str("*.out")))
                
            if len(out) < 1:
                print("Error: the directory", dirname, "has not *.out files!")
                sys.exit(1)

            else:
                xyz = get_xyz_gaussian(out)
                charge = get_charge(out)      #Charge
                mult = get_multiplicity(out)  #Multiplicity
                atom = xyz['atom']            #Atom
                x = xyz['x']                  #X coordinate in Angstrom
                y = xyz['y']                  #Y coordinate in Angstrom
                z = xyz['z']                  #Z coordinate in Angstrom  
                write_input(mem, keys, charge, mult, atom, x, y, z, dirname, out, unrestricted)
        else:
            print("too many arguments for this option!")
            sys.exit(0)
                
    if opt == '2':
        if len(args) == 6:
            out = sorted(get_output(str("*.out")))
                
            if len(out) < 1:
                print("Error: the directory", dirname, "has not *.out files!")
                sys.exit(1)

            else:   
                out = sorted(get_output(str("*.out")))
                irc = sort_by_rxcoord(out)
                rxcoord = irc['rxcoord']    #Reaction coordinate
                energy = irc['energy']      #Potential energy in Hartrees
                charge = irc['charge']      #Charge
                mult = irc['multiplicity']  #Multiplicity
                atom = irc['atom']          #Atom
                x = irc['x']                #X coordinate in Angstrom
                y = irc['y']                #Y coordinate in Angstrom
                z = irc['z']                #Z coordinate in Angstrom
                write_input_from_irc(mem, keys, charge, mult, atom, x, y, z, dirname, rxcoord, energy, unrestricted)
        else:
            print("too many arguments for this option!")
            sys.exit(0)

    if opt == '3':
        if len(args) >= 7:
            out = sorted(get_output(str("*.out")))
                
            if len(out) < 1:
                print("Error: the directory", dirname, "has not *.out files!")
                sys.exit(1)

            else:
                l = len(sys.argv[7])
                num_atoms = sys.argv[7][1:l-1].split(',')

                if len(num_atoms) < 1:
                    print("Error: atoms_list argument must be a list of numbers!")
                    sys.exit(1)

                else:
                    irc = sort_by_rxcoord(out)
                    rxcoord = irc['rxcoord']   #Reaction coordinate
                    energy = irc['energy']     #Potential energy in Hartrees
                    charge = irc['charge']     #Charge
                    mult = irc['multiplicity'] #Multiplicity
                    atom = irc['atom']         #Atom
                    x = irc['x']               #X coordinate in Angstrom
                    y = irc['y']               #Y coordinate in Angstrom
                    z = irc['z']               #Z coordinate in Angstrom
                    
                    natoms = len(atom[0])

                    list_num_atoms = [ ]
                    for i in num_atoms:
                        if int(i) <= 0:
                            print("Error: available first_atom argument is 1")
                            sys.exit(0)
                        if int(i) > (natoms):
                            print("Error: available last_atom argument is " + str(natoms))
                            sys.exit(0)
                        iatom = int(i)-1
                        list_num_atoms.append(int(iatom))

                        
                    all_atom_new = [ ]
                    all_x_new = [ ]
                    all_y_new = [ ]
                    all_z_new = [ ]
                    for (ia,ix,iy,iz) in zip(atom,x,y,z):
                        atom_new = [ ]
                        x_new = [ ]
                        y_new = [ ]
                        z_new = [ ]
                        for j in range(len(ia)):
                            for i in list_num_atoms:
                                if i == j:
                                    atom_new.append(ia[j])
                                    x_new.append(ix[j])
                                    y_new.append(iy[j])
                                    z_new.append(iz[j])
                                            
                        all_atom_new.append(atom_new)
                        all_x_new.append(x_new)
                        all_y_new.append(y_new)
                        all_z_new.append(z_new)
                                    
                    if len(args) == 8:
                        print("one argument is missing!")
                        sys.exit(0)
                                      
                    if len(args) == 9:
                        charge = int(sys.argv[8])
                        mult = int(sys.argv[9])

                    write_input_from_irc(mem, keys, charge, mult, all_atom_new, all_x_new, all_y_new, all_z_new, dirname, rxcoord, energy, unrestricted)

        elif len(args) < 7:
            print("some arguments are missing for this option!")
            sys.exit(0)
        elif len(args) > 9:
            print("too many arguments for this option!")
            sys.exit(0)

               
    if opt == '4':
        if len(args) >= 8:
            n_min = int(sys.argv[7]) - 1
            n_max = int(sys.argv[8]) - 1

            out = sorted(get_output(str("*.out")))
                
            if len(out) < 1:
                print("Error: the directory", dirname, "has not *.out files!")
                sys.exit(1)

            else:
                if n_min <= n_max:
                    if n_min < 0:
                        print("Error: available first_atom argument is 1")
                        sys.exit(0)
                        
                    else:
                        irc = sort_by_rxcoord(out)
                        rxcoord = irc['rxcoord']   #Reaction coordinate
                        energy = irc['energy']     #Potential energy in Hartrees
                        charge = irc['charge']     #Charge
                        mult = irc['multiplicity'] #Multiplicity
                        atom = irc['atom']         #Atom
                        x = irc['x']               #X coordinate in Angstrom
                        y = irc['y']               #Y coordinate in Angstrom
                        z = irc['z']               #Z coordinate in Angstrom

                        num_atoms = len(atom[0])
                        
                        if n_max > (num_atoms-1):
                            print("Error: available last_atom argument is " + str(num_atoms))
                            sys.exit(0)
                        
                        list_num_atoms = [ ]
                        list_num_atoms.extend(range(n_min, n_max))
                        list_num_atoms.append(n_max)

                        all_atom_new = [ ]
                        all_x_new = [ ]
                        all_y_new = [ ]
                        all_z_new = [ ]
                        for (ia,ix,iy,iz) in zip(atom,x,y,z):
                            atom_new = [ ]
                            x_new = [ ]
                            y_new = [ ]
                            z_new = [ ]
                            for j in range(len(ia)):
                                for i in list_num_atoms:
                                    if i == j:
                                        atom_new.append(ia[j])
                                        x_new.append(ix[j])
                                        y_new.append(iy[j])
                                        z_new.append(iz[j])
                                        
                            all_atom_new.append(atom_new)
                            all_x_new.append(x_new)
                            all_y_new.append(y_new)
                            all_z_new.append(z_new)
                                
                        if len(args) == 9:
                            print("one argument is missing!")
                            sys.exit(0)
                                  
                        if len(args) == 10:
                            charge = int(sys.argv[9])
                            mult = int(sys.argv[10])

                        write_input_from_irc(mem, keys, charge, mult, all_atom_new, all_x_new, all_y_new, all_z_new, dirname, rxcoord, energy, unrestricted)
                        
                else:
                    print("Error: last_atom argument must be greater-than first_atom argument")
                    sys.exit(0)

        elif len(args) < 8:
            print("some arguments are missing for this option!")
            sys.exit(0)

        elif len(args) > 10:
            print("too many arguments for this option!")
            sys.exit(0)

    if opt == '5':
        if len(args) == 8:
            xyz_files = sorted(get_output(str("*.xyz")))
            charge = int(sys.argv[7])
            mult = int(sys.argv[8])
                
            if len(xyz_files) < 1:
                print("Error: the directory", dirname, "has not *.xyz files!")
                sys.exit(1)

            else:
                xyz = read_xyz(xyz_files)
                atom = xyz['atom']            #Atom
                x = xyz['x']                  #X coordinate in Angstrom
                y = xyz['y']                  #Y coordinate in Angstrom
                z = xyz['z']                  #Z coordinate in Angstrom
                write_input_for_xyz(mem, keys, charge, mult, atom, x, y, z, dirname, xyz_files, unrestricted)
        else:
            print("too many arguments for this option!")
            sys.exit(0)

    if opt == '6':
        if len(args) == 6:
            out = sorted(get_output(str("*.out")))
                
            if len(out) < 1:
                print("Error: the directory", dirname, "has not *.out files!")
                sys.exit(1)

            else:   
                out = sorted(get_output(str("*.out")))
                scan_xyz = get_xyz_gaussian_scan(out)
                scan_out = get_scan_gaussian(out)
                charge = scan_out['charge'][0]      #Charge
                mult = scan_out['multiplicity'][0]  #Multiplicity
                atom = scan_xyz['atom'][0]             #Atom
                x = scan_xyz['x'][0]                   #X coordinate in Angstrom
                y = scan_xyz['y'][0]                   #Y coordinate in Angstrom
                z = scan_xyz['z'][0]                   #Z coordinate in Angstrom
                write_input_from_scan(mem, keys, charge, mult, atom, x, y, z, dirname, unrestricted)
        else:
            print("too many arguments for this option!")
            sys.exit(0)


    if opt == '7':
        if len(args) >= 7:
            out = sorted(get_output(str("*.out")))
                
            if len(out) < 1:
                print("Error: the directory", dirname, "has not *.out files!")
                sys.exit(1)

            else:
                l = len(sys.argv[7])
                num_atoms = sys.argv[7][1:l-1].split(',')

                if len(num_atoms) < 1:
                    print("Error: atoms_list argument must be a list of numbers!")
                    sys.exit(1)

                else:
                    out = sorted(get_output(str("*.out")))
                    scan_xyz = get_xyz_gaussian_scan(out)
                    scan_out = get_scan_gaussian(out)
                    charge = scan_out['charge'][0]      #Charge
                    mult = scan_out['multiplicity'][0]  #Multiplicity
                    atom = scan_xyz['atom'][0]             #Atom
                    x = scan_xyz['x'][0]                   #X coordinate in Angstrom
                    y = scan_xyz['y'][0]                   #Y coordinate in Angstrom
                    z = scan_xyz['z'][0]                   #Z coordinate in Angstrom
                    
                    natoms = len(atom[0])

                    list_num_atoms = [ ]
                    for i in num_atoms:
                        if int(i) <= 0:
                            print("Error: available first_atom argument is 1")
                            sys.exit(0)
                        if int(i) > (natoms):
                            print("Error: available last_atom argument is " + str(natoms))
                            sys.exit(0)
                        iatom = int(i)-1
                        list_num_atoms.append(int(iatom))

                        
                    all_atom_new = [ ]
                    all_x_new = [ ]
                    all_y_new = [ ]
                    all_z_new = [ ]
                    for (ia,ix,iy,iz) in zip(atom,x,y,z):
                        atom_new = [ ]
                        x_new = [ ]
                        y_new = [ ]
                        z_new = [ ]
                        for j in range(len(ia)):
                            for i in list_num_atoms:
                                if i == j:
                                    atom_new.append(ia[j])
                                    x_new.append(ix[j])
                                    y_new.append(iy[j])
                                    z_new.append(iz[j])
                                            
                        all_atom_new.append(atom_new)
                        all_x_new.append(x_new)
                        all_y_new.append(y_new)
                        all_z_new.append(z_new)
                                    
                    if len(args) == 8:
                        print("one argument is missing!")
                        sys.exit(0)
                                      
                    if len(args) == 9:
                        charge = int(sys.argv[8])
                        mult = int(sys.argv[9])

                    write_input_from_scan(mem, keys, charge, mult, all_atom_new, all_x_new, all_y_new, all_z_new, dirname, unrestricted)

        elif len(args) < 7:
            print("some arguments are missing for this option!")
            sys.exit(0)
        elif len(args) > 9:
            print("too many arguments for this option!")
            sys.exit(0)

    if opt == '8':
        if len(args) >= 8:
            n_min = int(sys.argv[7]) - 1
            n_max = int(sys.argv[8]) - 1

            out = sorted(get_output(str("*.out")))
                
            if len(out) < 1:
                print("Error: the directory", dirname, "has not *.out files!")
                sys.exit(1)

            else:
                if n_min <= n_max:
                    if n_min < 0:
                        print("Error: available first_atom argument is 1")
                        sys.exit(0)
                        
                    else:
                        out = sorted(get_output(str("*.out")))
                        scan_xyz = get_xyz_gaussian_scan(out)
                        scan_out = get_scan_gaussian(out)
                        charge = scan_out['charge'][0]      #Charge
                        mult = scan_out['multiplicity'][0]  #Multiplicity
                        atom = scan_xyz['atom'][0]             #Atom
                        x = scan_xyz['x'][0]                   #X coordinate in Angstrom
                        y = scan_xyz['y'][0]                   #Y coordinate in Angstrom
                        z = scan_xyz['z'][0]                   #Z coordinate in Angstrom

                        num_atoms = len(atom[0])
                        
                        if n_max > (num_atoms-1):
                            print("Error: available last_atom argument is " + str(num_atoms))
                            sys.exit(0)
                        
                        list_num_atoms = [ ]
                        list_num_atoms.extend(range(n_min, n_max))
                        list_num_atoms.append(n_max)

                        all_atom_new = [ ]
                        all_x_new = [ ]
                        all_y_new = [ ]
                        all_z_new = [ ]
                        for (ia,ix,iy,iz) in zip(atom,x,y,z):
                            atom_new = [ ]
                            x_new = [ ]
                            y_new = [ ]
                            z_new = [ ]
                            for j in range(len(ia)):
                                for i in list_num_atoms:
                                    if i == j:
                                        atom_new.append(ia[j])
                                        x_new.append(ix[j])
                                        y_new.append(iy[j])
                                        z_new.append(iz[j])
                                        
                            all_atom_new.append(atom_new)
                            all_x_new.append(x_new)
                            all_y_new.append(y_new)
                            all_z_new.append(z_new)
                                
                        if len(args) == 9:
                            print("one argument is missing!")
                            sys.exit(0)
                                  
                        if len(args) == 10:
                            charge = int(sys.argv[9])
                            mult = int(sys.argv[10])

                        write_input_from_scan(mem, keys, charge, mult, all_atom_new, all_x_new, all_y_new, all_z_new, dirname, unrestricted)
                        
                else:
                    print("Error: last_atom argument must be greater-than first_atom argument")
                    sys.exit(0)

        elif len(args) < 8:
            print("some arguments are missing for this option!")
            sys.exit(0)

        elif len(args) > 10:
            print("too many arguments for this option!")
            sys.exit(0)
                
    if opt == '9':
        if len(args) == 6:
            out = sorted(get_output(str("*.out")))
                
            if len(out) < 1:
                print("Error: the directory", dirname, "has not *.out files!")
                sys.exit(1)

            else:   
                out = sorted(get_output(str("*.out")))
                irc = sort_by_rxcoord(out)
                rxcoord = irc['rxcoord']    #Reaction coordinate
                energy = irc['energy']      #Potential energy in Hartrees
                charge = irc['charge']      #Charge
                mult = irc['multiplicity']  #Multiplicity
                atom = irc['atom']          #Atom
                x = irc['x']                #X coordinate in Angstrom
                y = irc['y']                #Y coordinate in Angstrom
                z = irc['z']                #Z coordinate in Angstrom
                write_inputchk_from_irc(mem, keys, charge, mult, atom, x, y, z, dirname, rxcoord, energy, unrestricted)
        else:
            print("too many arguments for this option!")
            sys.exit(0)   
        
        
#----------------------------------------------------------

a = write_input()
