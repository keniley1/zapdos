#Find "C", going line by line, if C is found, check if character next
#to it is an integer - if it is, it is a reaction, read the rate coeff
#data that is 2 lines below it into a file and name that file the C# reaction
#designation. Add that C# designation to a list, next time it searches for 
#C, check the C# against the list, if it is there then skip and continue 
#searching the file. Use file.seek(0) to reset to beginning 
#(Or do file.seek(last_reac) where last_reac is the line number of the last
#reaction it read, and start from there. This could eliminate need to check if
#it has done the reaction already
#Repeat until end of bolsig datafile is reached.
#Function to create a datafile with the reaction designation as its name
import numpy as np

def create_k_file(reaction_name):
    filename = "%s.txt" % reaction_name
    f = open(filename,"w+")
    f.close()
    return filename

def read_bolsig(bolsigdatafile, XSdatafile, runs, gas_density, coefficient_type="rate", transport_factor=1.0, rate_factor=1.0): #Takes name of bolsig data file as input
    '''
    transport_factor : factor which multiplies the transport coefficients to scale the units.
    rate_factor      : factor which multiplies the rate coefficients to scale the units.
    '''
    with open(XSdatafile,"r") as fXS:
        rxns = []
        for line in fXS:
            if line.find("PROCESS:") != -1:
                rxns.append(line[11:line.find(",")])

    with open(bolsigdatafile,"r") as fbolsig:
        lines = fbolsig.readlines()
        reac_lines = [] #list of the line numbers with reactions, to be updated with each reaction read.
        kfile_names = [] #names listed as C1 ... CN
        #kfile_names = rxns #list of the file names (C1.txt...CN.txt)
        e_mob_dif_data = []
        mobility = []
        diffusivity = []
        Efield = []
        reac_num = 0
        for i in range(len(lines)):#Find lines with the rxns, make datafiles for each
            line = lines[i]
            isreac = line.find("C") #gives -1 if "C" is not in line
            #reac_num = 0 #number of reactions found in the bolsig file
            if isreac != -1: # If C is in line
                #Check if C is followed by a digit. if so, this is a reaction, read its data into a file
                if (line[isreac+1].isdigit() and line[isreac-1] != "("):
                    #print(line.replace(" ", "_"))
                    temp_name = line.rstrip()
                    temp_name = temp_name.replace(" ", "_")
                    temp_name = temp_name.replace("__", "_")
                    temp_name = temp_name.replace("__", "_")
                    temp_name = temp_name.replace("(", "")
                    temp_name = temp_name.replace(")", "")
                    kfile_names.append(temp_name.replace("\n", ""))
#                    reacname = line.split()[0]
#                    kfilename = create_k_file(reacname)
#                    kfilename = create_k_file(kfile_names[reac_num])
                    reac_num += 1
                    reac_lines.append(i)
#                    kfile_names.append(kfilename)
            if line.find("Mobility") != -1:
                    for k in range(i+1,i+runs+1):
                        mobility.append(lines[k].split()) # list of pairs of energy and mobility data
            if line.find("Diffusion") != -1:
                num = 0
                mu = 0
                diff = 0
                with open("electron_mobility_diffusion.txt","w+") as f:
                    for k in range(i+1,i+runs+1):
                        diffusivity.append(lines[k].split())

                        mu = str(float(mobility[num][1]) * transport_factor / gas_density)
                        diff = str(float(diffusivity[num][1]) * transport_factor / gas_density)
                        #f.write("%s\t%s\t%s\n" % (diffusivity[num][0], mobility[num][1], diffusivity[num][1]))
                        f.write("%s\t%s\t%s\n" % (diffusivity[num][0], mu, diff))
                        num += 1
            if line.find("Electric") != -1:
                num = 0
                #with open("reduced_field.txt", "w+") as f:
                for k in range(i+1, i+runs+1):
                    Efield.append(lines[k].split())

#Now reac_lines has the line numbers of the first header for each reaction dataset
#data columns begin 2 lines down from these headers, and are 300 lines long each, with 2 
#lines of whitespace between end of dataset and beginning of next dataset
        for j in range(len(reac_lines)):
            #f = open(kfile_names[j],"w+")
            x_data = np.zeros(shape=(runs,))
            y_data = np.zeros(shape=(runs,))
            cc = 0
            tstr = "C12"
            for i in range(reac_lines[j]+2,reac_lines[j]+runs+2):
                tempx, tempy = lines[i].split()
                x_data[cc] = float(tempx)
                y_data[cc] = float(tempy)
                cc += 1

                #f.write(lines[i])
            if (coefficient_type == "townsend"):
                with open(kfile_names[j],'w+') as write_file:
                    for i in range(runs):
                        mu = float(mobility[i][1])
                        E_N = float(Efield[i][1])*1e-21
                        write_file.write('{0:.4e} {1:.4e}\n'.format(x_data[i], y_data[i]*6.022e23 / (mu * E_N ) * transport_factor))
            else:
                with open(kfile_names[j],'w+') as write_file:
                    for i in range(runs):
                        #write_file.write('{0:.4e} {1:.4e}\n'.format(x_data[i], y_data[i]*1e6))
                        write_file.write('{0:.4e} {1:.4e}\n'.format(x_data[i], y_data[i]*6.022e23 * rate_factor))
                #f.close()
            
    return None

#read_bolsig("bolsigdb_out.dat", "argon_cs_data.dat", runs=100, gas_density=2.445692e25, coefficient_type="townsend", transport_factor=1, rate_factor=1)
read_bolsig("bolsigdb_out.dat", "total_cross_sections_paschen.txt", runs=100, gas_density=2.445692e25, coefficient_type="townsend", transport_factor=1, rate_factor=1)
