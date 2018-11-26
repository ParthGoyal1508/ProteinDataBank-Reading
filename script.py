import numpy as np
import math


def get_dihedral_angle(a, b, c, d):
    v1 = (b[0]-a[0], b[1]-a[1], b[2]-a[2])
    v2 = (c[0]-b[0], c[1]-b[1], c[2]-b[2])
    v3 = (d[0]-c[0], d[1]-c[1], d[2]-c[2])
    normal1 = np.cross(v1, v2)
    normal2 = np.cross(v2, v3)
    # unit vector along normal n1
    normal1 = (np.cross(v1, v2)/np.linalg.norm(normal1))
    # unit vector along normal n2
    normal2 = (np.cross(v2, v3)/np.linalg.norm(normal2))
    unit_v2 = (v2/np.linalg.norm(v2))  # unit vector along normal v2
    m1 = (np.cross(normal1, unit_v2))
    x = normal1[0]*normal2[0] + normal1[1]*normal2[1] + normal1[2]*normal2[2]
    y = m1[0]*normal2[0] + m1[1]*normal2[1] + m1[2]*normal2[2]
    return(round(math.degrees(math.atan2(y, x)), 3))


def main():
    inp = input()
    fp = open(inp, "r")
    out = open(inp.split('.')[0]+'_output.txt', 'w')

    array = fp.readlines()

    chain = {}
    aminoacid = {}
    ligand = {}
    li = []
    angles = []
    temp1 = {}
    temp2 = {}
    temp1['N'] = temp2['N'] = (float('nan'), float('nan'), float('nan'))
    temp1['CA'] = temp2['CA'] = (float('nan'), float('nan'), float('nan'))
    temp1['C'] = temp2['C'] = (float('nan'), float('nan'), float('nan'))

    angle_index = 0
    angle_tuple = (0, 0, 0)
    cur_chain = None
    sumoflength = 0
    unk = 0
    name = ""
    chainid = ""

    i = -1

    while i < len(array)-1:
        i += 1
        line = array[i]
        if line.split(' ',)[0] == 'TITLE':
            line1 = line.strip()  # remove last space
            line1 = ' '.join(line1.split())
            name += line1.split(' ', 1)[1]
        if line.split(' ',)[0] == 'SEQRES':
            line = line.strip()
            var = line.split()
            x = var[2]
            if x not in chain.keys():
                chain[x] = int(var[3])
                sumoflength += chain[x]
            y = var[4:]
            for j in y:
                if j not in aminoacid:
                    aminoacid[j] = 0
                aminoacid[j] += 1
        if line.split(' ',)[0] == 'HET':
            line1 = line.strip()
            line1 = " ".join(line1.split())
            line1 = line1.split(" ")
            lig = line1[1]
            if lig not in ligand.keys():
                ligand[lig] = 0
            ligand[lig] += 1
        if line.split(' ',)[0] == 'ATOM':
            line1 = line.strip()
            line1 = " ".join(line1.split())
            line1 = line1.split(" ")
            if line1[4] == cur_chain:
                if angle_index == 1:
                    temp1['N'] = temp2['N']
                    temp1['CA'] = temp2['CA']
                    temp1['C'] = temp2['C']
                    while True:
                        if line1[0] != 'ATOM' or line1[4] != cur_chain:
                            angles.append(('psi', float('nan')))
                            angles.append(('omega', float('nan')))
                            angle_index = (angle_index + 2) % 3
                            break
                        if line1[2] == 'N':
                            temp2['N'] = (float(line1[6]), float(
                                line1[7]), float(line1[8]))
                            psi = get_dihedral_angle(
                                temp1['N'], temp1['CA'], temp1['C'], temp2['N'])
                            angle_tuple = ('psi', psi)
                            angles.append(angle_tuple)
                            angle_index = (angle_index + 1) % 3
                            break
                        i += 1
                        line = array[i]
                        line1 = line.strip()
                        line1 = " ".join(line1.split())
                        line1 = line1.split(" ")
                    continue

                if angle_index == 2:
                    temp2['CA'] = (float(line1[6]), float(
                        line1[7]), float(line1[8]))
                    omega = get_dihedral_angle(
                        temp1['CA'], temp1['C'], temp2['N'], temp2['CA'])
                    angle_tuple = ('omega', omega)
                    angles.append(angle_tuple)
                    angle_index = (angle_index + 1) % 3
                    continue

                if angle_index == 0:
                    temp2['C'] = (float(line1[6]), float(
                        line1[7]), float(line1[8]))
                    phi = get_dihedral_angle(
                        temp1['C'], temp2['N'], temp2['CA'], temp2['C'])
                    angle_tuple = ('phi', phi)
                    angles.append(angle_tuple)
                    angle_index = (angle_index + 1) % 3
                    continue

            else:
                cur_chain = line1[4]
                temp2['N'] = (float(line1[6]), float(
                    line1[7]), float(line1[8]))
                i += 1
                line = array[i]
                if line.split(' ',)[0] != 'ATOM':
                    i += 1
                    line = array[i]
                line1 = line.strip()
                line1 = " ".join(line1.split())
                line1 = line1.split(" ")
                temp2['CA'] = (float(line1[6]), float(
                    line1[7]), float(line1[8]))
                i += 1
                line = array[i]
                if line.split(' ',)[0] != 'ATOM':
                    i += 1
                    line = array[i]
                line1 = line.strip()
                line1 = " ".join(line1.split())
                line1 = line1.split(" ")
                temp2['C'] = (float(line1[6]), float(
                    line1[7]), float(line1[8]))

                angle_tuple = ('phi', float('nan'))
                angles.append(angle_tuple)
                angle_index = (angle_index + 1) % 3

    out.write(name+'\n')
    out.write("LENGTH"+'\t'+str(sumoflength)+'\n')
    out.write("CHAINS"+'\t'+str(len(chain.keys())) +
              '\t'+','.join(sorted(chain.keys()))+'\n')
    for k in sorted(aminoacid.keys()):
        if k == 'UNK':
            unk += 1
        else:
            out.write(k+'\t'+str(float(aminoacid[k])/sumoflength)+'\n')
    if unk == 0:
        out.write("UNKNOWN"+'\t'+"0"+'\n')
    if unk != 0:
        out.write("UNKNOWN"+'\t'+str(aminoacid["UNK"])+'\n')
    if len(ligand.keys()) != 0:
        out.write("LIGANDS"+'\t'+','.join(sorted(ligand.keys()))+'\n')
    i = 0
    for k in sorted(chain.keys()):
        out.write('CHAIN-' + str(k)+'\n')
        while i < len(angles)-1:
            phi = angles[i][1]
            psi = angles[i+1][1]
            omega = angles[i+2][1]
            if math.isnan(phi):
                phi = 'NA'
            if math.isnan(psi):
                psi = 'NA'
            if math.isnan(omega):
                omega = 'NA'
            out.write(str(phi)+'\t'+str(psi)+'\t'+str(omega)+'\n')
            i += 3
            if psi == 'NA' and omega == 'NA':
                break
    fp.close()


if __name__ == '__main__':
    main()
