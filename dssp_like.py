#!/usr/bin/python

#-------PROJET  Assignation des structures secondaire de protéines ------- #
#-------------------- TIEO SONIA - M2BI -----------------------------------#




#*********************IMPORTATIONS *********************#
import sys , getopt
import math
import numpy as np
from itertools import groupby
from operator import itemgetter



#*********************GLOBAL VARIABLES*********************#

#List of targeted atoms ( electrostatic interaction )
interac_atoms = ["N", "C", "O" , "H" , "CA"]
# q1, q2 : partial charges
q1 = 0.42
q2 = 0.20
#dimensional factor
f = 332
# Energy cutoff for Hbond
Emin = -0.5



#*********************FONCTIONS*********************#

def extract_atom_coord(line, dico_atom) :
    """Extracts cartesian coordinates, residue number and name of an atom

    Parameters
    ----------
    line : str
        line of this atom in PDB
    dico_atom : dict
        of an atom in a PDB

    Returns
    -------
    dict
        dict of this atom with informations

    """
    atom = dico_atom
    atom['res_num'] = int(line[22:26])
    atom['res_name'] = line[17:20].strip()
    atom['x'] = float(line[30:38])
    atom['y'] = float(line[38:46])
    atom['z'] = float(line[46:54])
    return(atom)


def extract_first_resid(pdbfile):
    """Extracts number of first residue

    Parameters
    ----------
    pdbfile : file with PDB format

    Returns
    -------
    int
        number of first residue

    """
    with open( pdbfile, "r") as pdb_file:
        for line in pdb_file:
            atom_name = line[12:16].strip()
            if (line.startswith("ATOM") and atom_name in interac_atoms ):
                res_first = int(line[22:26])
                break
    return(res_first)


def parse_pdb(pdbfile):
    """Short summary.

    Parameters
    ----------
    pdbfile : file with PDB format

    Returns
    -------
    dict
        dict of Dict of all residues - key : residue number , item : dict
        for each residue - key : atom name , item : coordinates

    """

    with open( pdbfile, "r") as pdb_file:
        res = extract_first_resid(pdbfile) - 1
        res_list = []
        resid = {}
        for line in pdb_file:
            atom_name = line[12:16].strip()
            if (line.startswith("ATOM") and atom_name in interac_atoms ):
                res_num = int(line[22:26])
                if res_num in res_list:
                    if (atom_name not in resid[res]) :
                        atom = {}
                        resid[res][atom_name] = extract_atom_coord(line, atom)
                    else:
                        resid[res][atom_name] = extract_atom_coord(line, atom)
                else :
                    res_list.append(res_num)
                    res +=1
                    resid[res] = {}
                    atom = {}
                    resid[res][atom_name] = extract_atom_coord(line, atom)

    return(resid)

def dist2atoms(dico_res, numresA, numresB, atomA, atomB ):
    """Extracts distance between 2 atoms from coordinates

    Parameters
    ----------
    dico_res : dict
        dict with all residues
    numresA : int
        residue A number
    numresB : int
        residue B number
    atomA : str
        atom A name like 'C' , 'N' ...
    atomB : str
        atom B name like 'C' , 'N' ....

    Returns
    -------
    int
        distance angstrom units

    """
    xA = dico_res[numresA][atomA]['x']
    yA = dico_res[numresA][atomA]['y']
    zA = dico_res[numresA][atomA]['z']
    xB = dico_res[numresB][atomB]['x']
    yB = dico_res[numresB][atomB]['y']
    zB = dico_res[numresB][atomB]['z']
    rayon = math.sqrt(( xB-xA)**2 + ( yB-yA)**2 + ( zB-zA)**2 )
    return(rayon)



def is_Hbond(dico_res, numresA, numresB):
    """Hbond (A,B) = Electrostatic interaction energy between 2 residue A, B  < 0.5kcal/mole

    Parameters
    ----------
    dico_res : dict
        dict with all residues
    numresA : int
        residue A number
    numresB : int
        residue B number

    Returns
    -------
    bool
        return True if HBond(A,B)

    """
    rON = dist2atoms(dico_res , numresA , numresB , 'O' , 'N' )
    rOH = dist2atoms(dico_res , numresA , numresB , 'O' , 'H' )
    rCN = dist2atoms(dico_res , numresA , numresB , 'C' , 'N' )
    rCH = dist2atoms(dico_res , numresA , numresB , 'C' , 'H' )
    E_elec = q1 * q2 * f * (1/rON + 1/rCH - 1/rOH - 1/rCN)
    if (E_elec < Emin ) :
        return ( True)
    return(False)


def list_nturn(dico_res, nturn ):
    """

    Parameters
    ----------
    dico_res : dict
        dict with all residues
    nturn : int
        nturn = (3,4 or 5)

    Returns
    -------
    list of list
        list[0] : all the residues involved in n-list_turn
        list[1] : list of tuples of 2 bonding residue

    """
    list_turn = []
    list_partner = []
    res_first = list ( dico_res.keys())[0]
    size = len(list ( dico_res.keys() ))
    for res in range ( res_first , ( size + res_first -1 ) - nturn   ) :
        if ( 'H' in dico_res[res+nturn ].keys()) :
            if is_Hbond( dico_res , res, res+nturn ) :
                list_turn.append(res)
                list_turn.append(res+nturn)
                list_partner.append((res , res+nturn))
        if len(list_partner) > 0 :
            first,snd = zip(*list_partner)
    return( list(set(list_turn) ), list_partner)



def successives_nb(liste):
    """ Finds succssives numbers in a given list

    Parameters
    ----------
    liste : list
        list of int

    Returns
    -------
    list
        list on int

    """
    succ_nb = []
    for k,g in groupby(enumerate(sorted(liste)),lambda x:x[0]-x[1]):
        group = (map(itemgetter(1),g))
        group = list(map(int,group))
        if len(group) > 1 :
            succ_nb = succ_nb + group
    return(succ_nb)


def helix_alpha(list_4turn):
    """ Finds residues in alpha helix from residues in 4-turn list of list

    Parameters
    ----------
    list_4turn : list of list
        list_4turn[0] : : all the residues involved in 4-list_turn

    Returns
    -------
    list
        list of residues in alpha helix

    """
    list_ha = successives_nb (list_4turn[0])
    return (list_ha)


def helix_310 (list_3turn) :
    """ Finds residues in 3,10 helix from residues in 3-turn list of list

    Parameters
    ----------
    list_3turn : list of list
        list_3turn[0] : : all the residues involved in 3-list_turn

    Returns
    -------
    list
        list of residues in 3,10 helix

    """
    list_h310 = successives_nb (list_3turn[0])
    return (list_h310)


def helix_pi (list_5turn):
    """ Finds residues in pi helix from residues in 5-turn list of list

    Parameters
    ----------
    list_5turn : list of list
        list_5turn[0] : : all the residues involved in 5-list_turn

    Returns
    -------
    list
        list of residues in pi helix

    """
    list_hpi = successives_nb (list_5turn[0])
    return (list_hpi)


def parallel_bridge(dico_res):
    """Finds residues in parallel bridge

    Parameters
    ----------
    dico_res : dict
        dict with all residues

    Returns
    -------
    list of list
        list of all residues i , list of all i bridge parteners j

    """
    res_first = list ( dico_res.keys())[0]
    size = len(list ( dico_res.keys()))
    list_i =[]
    list_j = []
    for i in range( res_first + 1 , ( size + res_first - 1 ) -1 ):
        for j in range( i+3 , ( size + res_first - 1 ) -1  ) :
            if ( 'H' in dico_res[j].keys()) and ( 'H' in dico_res[i+1].keys()):
                cas1 = is_Hbond(dico_res, i-1, j) and is_Hbond(dico_res, j, i+1)
            else :
                cas1=False
            if ( 'H' in dico_res[i].keys()) and ( 'H' in dico_res[j+1].keys()):
                cas2 = is_Hbond(dico_res, j-1, i) and is_Hbond(dico_res, i, j+1)
            else :
                cas2=False
            if cas1==True or cas2==True :
                list_i.append(i)
                list_j.append(j)
    return([list_i , list_j])


def antiparallel_bridge(dico_res):
    """Finds residues in antiparallel bridge

    Parameters
    ----------
    dico_res : dict
        dict with all residues

    Returns
    -------
    list of list
        list of all residues i , list of all i bridge parteners j

    """
    res_first = list ( dico_res.keys())[0]
    size = len(list ( dico_res.keys()))
    list_i =[]
    list_j = []
    indexes=[]
    for i in range( res_first + 1 , ( size + res_first - 1 )  -1 ):
        for j in range(  i + 3 , (size + res_first - 1 ) - 1   ) :
            if ( 'H' in dico_res[j].keys()) and ( 'H' in dico_res[i].keys()):
                cas1 = is_Hbond(dico_res, i, j) and is_Hbond(dico_res, j, i)
            else :
                cas1 =False
            if ( 'H' in dico_res[j-1].keys()) and ( 'H' in dico_res[i+1].keys() ) and ('H' in dico_res[j+1].keys()) and ('H' in dico_res[i-1].keys()):
                cas2 = is_Hbond(dico_res, i-1, j+1) and is_Hbond(dico_res, j-1, i+1)
            else:
                cas2=False
            if cas1==True or cas2==True :
                list_i.append(i)
                list_j.append(j)
    return([list_i , list_j ])


def partner_bridge( list_ij ) :
    """Finds bridge partner(s) i,j

    Parameters
    ----------
    list_ij : list of list
        list of all residues i , list of all i bridge parteners j

    Returns
    -------
    list of tuples
        list of partners bridge (i,j)

    """
    partner = []
    list_i = list_ij[0]
    list_j = list_ij[1]
    for i in range(len(list_i)) :
        partner.append( (list_i[i] , list_j[i]))
    for j in range(len(list_j)) :
        partner.append(( list_j[j] , list_i[j]))
    return(partner)


def vect2atoms(dico_res, numresA, numresB, atomA, atomB) :
    """ Vector AB from coordinates of atom A and atom B

    Parameters
    ----------
    dico_res : dict
        dict with all residues
    numresA : int
        residue A number
    numresB : int
        residue B number
    atomA : str
        atom A name like 'C' , 'N' ...
    atomB : str
        atom B name like 'C' , 'N' ....

    Returns
    -------
    tuple
        vectAB( xB-xA , yB-yA ,zB-zA)

    """
    xA = dico_res[numresA][atomA]['x']
    yA = dico_res[numresA][atomA]['y']
    zA = dico_res[numresA][atomA]['z']
    xB = dico_res[numresB][atomB]['x']
    yB = dico_res[numresB][atomB]['y']
    zB = dico_res[numresB][atomB]['z']
    vectAB = ( xB-xA , yB-yA ,zB-zA)
    return(vectAB)



def unit_vector(vector):
    """Returns the norm of vector from the coordinates of vector


    Parameters
    ----------
    vector : tuple

    Returns
    -------
    numpy.ndarray
        norm of vector

    """
    return (vector / np.linalg.norm(vector))

def angle_between(v1, v2):
    """Returns the angle in degrees between vectors v1 and v2


    Parameters
    ----------
    v1 : tuple
    v2 : tuple

    Returns
    -------
    float
        angle in degrees

    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return (math.degrees(np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))) )


def bend( dico_res ) :
    """Finds residues involved in bends

    Parameters
    ----------
    dico_res : dict
        dict with all residues

    Returns
    -------
    list
        list on residue numbers

    """
    list_bend=[]
    atomA = atomB = 'CA'
    res_first = list ( dico_res.keys())[0]
    size = len(list ( dico_res.keys()))
    for res in range( res_first + 2 , ( size + res_first - 1 )  - 2 ):
        v1 = vect2atoms(dico_res , res - 2 , res ,'CA' , 'CA' )
        v2 = vect2atoms(dico_res , res , res + 2 ,'CA' , 'CA' )
        if angle_between(v1 , v2) > 70 :
            list_bend.append(res)
    return(list_bend)

def dict_resid_structure ( dico_res ) :
    """Creates "empty" dict with structural labels for all residues

    Parameters
    ----------
    dico_res : dict
        dict with all residues

    Returns
    -------
    dict
        empty with structural labels for all residues

    """
    res_struct = {}
    for res in dico_res :
        struct ={}
        res_struct[res] = struct
        struct['res_num'] =  dico_res[res]['N']['res_num']
        struct['res_name'] = dico_res[res]['N']['res_name']
        struct['type'] = ""
        struct['BP1'] = 0
        struct['BP2'] = 0
        struct['3T'] = ""
        struct['4T'] = ""
        struct['5T'] = ""
    return(res_struct)


def res_bend ( dico_res ,  dico_res_struct , list_bend) :
    """Refills dict with structural informations for all residues

    Parameters
    ----------
    dico_res : dict
        dict with all residues
    dico_res_struct : dict
        dict with all residues structures informations
    list_bend : list
        residues involved in bend

    Returns
    -------
    dict
        dico_res_struct

    """
    for res in dico_res:
        if res in list_bend  :
            dico_res_struct[res]['type'] = 'S'
    return( dico_res_struct)

def res_turnT( dico_res , dico_res_struct , list_4turn , list_3turn, list_5turn ) :
    """Refills dict with structural informations for all residues

    Parameters
    ----------
    dico_res : dict
        dict with all residues
    dico_res_struct : dict
        dict with all residues structures informations
    list_4turn : list of list
        list_4turn[0] : : all the residues involved in 4-list_turn
    list_3turn : list of list
        list_3turn[0] : : all the residues involved in 3-list_turn
    list_5turn : list of list
        list_5turn[0] : : all the residues involved in 4-list_turn

    Returns
    -------
    dict
        dico_res_struct
    """
    for res in dico_res :
        if res in (set(list_4turn[0])) or res in (set(list_5turn[0]) )  or res in (set(list_3turn[0])) :
            dico_res_struct[res]['type'] = 'T'
        if res in [x[0] for x in list_3turn[1]] :
            dico_res_struct[res]["3T"]= '>'
            dico_res_struct[res + 1]["3T"]= '3'
            dico_res_struct[res + 2]["3T"]= '3'
        if res in [x[1] for x in list_3turn[1]] :
            dico_res_struct[res]['3T'] = '<'
        if res in [x[0] for x in list_4turn[1]] :
            dico_res_struct[res]["4T"]= '>'
            dico_res_struct[res + 1]["4T"]= '4'
            dico_res_struct[res + 2]["4T"]= '4'
            dico_res_struct[res + 3]["4T"]= '4'
        if res in [x[1] for x in list_4turn[1]] :
            dico_res_struct[res]['4T'] = '<'
        if res in [x[0] for x in list_5turn[1]] :
            dico_res_struct[res]["5T"]= '>'
            dico_res_struct[res + 1]["5T"]= '5'
            dico_res_struct[res + 2]["5T"]= '5'
            dico_res_struct[res + 3]["5T"]= '5'
            dico_res_struct[res + 4]["5T"]= '5'
        if res in [x[1] for x in list_5turn[1]] :
            dico_res_struct[res]['5T'] = '<'
    return ( dico_res_struct)



def res_betasheetE ( dico_res ,  dico_res_struct ,  par , antipar ) :
    """Refills dict with structural informations for all residues

    Parameters
    ----------
    dico_res : dict
        dict with all residues
    dico_res_struct : dict
        dict with all residues structures informations
    par : list of list
        list of all residues i , list of all i bridge parteners j
    antipar : list of list
        list of all residues i , list of all i bridge parteners j

    Returns
    -------
    dict
        dico_res_struct

    """
    for res in dico_res :
        if res in par[0] or res in par[1] :
            dico_res_struct[res]['type'] = 'E'
            partner = partner_bridge(par)
            list_partner = [item for item in partner if item[0] == res]
            dico_res_struct[res]['BP1']  = list_partner[0][1]
            if len(list_partner) == 2 :
                dico_res_struct[res]['BP2']  = list_partner[1][1]
        if res in antipar[0] or res in antipar[1] :
            dico_res_struct[res]['type'] = 'E'
            partner = partner_bridge(antipar)
            list_partner = [item for item in partner if item[0] == res]
            dico_res_struct[res]['BP1']  = list_partner[0][1]
            if len(list_partner) == 2 :
                dico_res_struct[res]['BP2']  = list_partner[1][1]
    return( dico_res_struct)



def res_helixpi ( dico_res ,  dico_res_struct ,  h_pi) :
    """Refills dict with structural informations for all residues

    Parameters
    ----------
    dico_res : dict
        dict with all residues
    dico_res_struct : dict
        dict with all residues structures informations
    h_pi : list
        Description of parameter `h_pi`.

    Returns
    -------
    dict
        dico_res_struct

    """
    for res in dico_res :
        if res in h_pi :
            dico_res_struct[res]['type'] = 'I'
    return(dico_res_struct)


def res_helixH ( dico_res ,  dico_res_struct ,  h_alpha) :
    """Refills dict with structural informations for all residues

    Parameters
    ----------
    dico_res : dict
        dict with all residues
    dico_res_struct : dict
        dict with all residues structures informations
    h_alpha : list
        list of residues in alpha helix

    Returns
    -------
    dict
        dico_res_struct

    """
    for res in dico_res :
        if res in h_alpha :
            dico_res_struct[res]['type'] = 'H'
    return(dico_res_struct)



def res_helix310 ( dico_res ,  dico_res_struct ,  h_310) :
    """Refills dict with structural informations for all residues

    Parameters
    ----------
    dico_res : dict
        dict with all residues
    dico_res_struct : dict
        dict with all residues structures informations
    h_310 : list
        list of residues in 3,10 helix

    Returns
    -------
    dict
        dico_res_struct

    """
    for res in dico_res :
        if res in h_310 :
            dico_res_struct[res]['type'] = 'G'
    return(dico_res_struct)


def affichage( dico_res_struct , outputfile ) :
    """Write structures assignation for each residue

    Parameters
    ----------
    dico_res_struct : dict
        dico_res_struct


    """
    with open(outputfile, 'w') as fillout:
        fillout.write("{:>3s}{:>4s}{:>12s}{:>4s}{:>4s}\n".format( "RES", "AA" , "STRUCTURE" , "BP1", "BP2"))
        for res in dico_res_struct :
            r = dico_res_struct[res]
            fillout.write("{:>3d}{:>4s}{:>3s}{:>3s}{:>3s}{:>3s}{:>4d}{:>4d}\n".format(res, r['res_name'] , r['type'], r['3T'] , r['4T'] , r['5T'] , r['BP1'] , r['BP2']))




#********************* MAIN *********************#


def main(argv):
   pdbfile = ''
   outputfile = ''
   try:
      opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
   except getopt.GetoptError:
      print ('dssp_like.py -i <pdbfile> -o <outputfile>')
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print ('dssp_like.py -i <pdbfile> -o <outputfile>')
         sys.exit()
      elif opt in ("-i", "--ifile"):
         pdbfile = arg
      elif opt in ("-o", "--ofile"):
         outputfile = arg
   print ('Input file is', pdbfile )
   print ('Output file is', outputfile)
   dico_res = parse_pdb(pdbfile)
   list_4turn = (list_nturn( dico_res , 4 ) )
   list_5turn = (list_nturn( dico_res , 5 ) )
   list_3turn = (list_nturn( dico_res , 3 ) )
   h_310 = helix_310(list_3turn)
   h_alpha = helix_alpha(list_4turn)
   h_pi = (helix_pi(list_5turn))
   par = (parallel_bridge (dico_res))
   antipar = (antiparallel_bridge (dico_res))
   list_bend = bend( dico_res )
   dico_res_struct = dict_resid_structure ( dico_res  )
   dico_res_struct = res_turnT( dico_res , dico_res_struct , list_4turn , list_3turn, list_5turn )
   dico_res_struct = res_bend( dico_res , dico_res_struct , list_bend)
   dico_res_struct = res_betasheetE ( dico_res ,  dico_res_struct ,  par , antipar )
   dico_res_struct = res_helix310( dico_res ,  dico_res_struct , h_310)
   dico_res_struct = res_helixH( dico_res ,  dico_res_struct , h_alpha)
   dico_res_struct = res_helixpi( dico_res ,  dico_res_struct , h_pi)
   affichage( dico_res_struct , outputfile )
if __name__ == "__main__":
   main(sys.argv[1:])

#pdbfile = str( sys.argv[1] )
#outputfile = str( sys.argv[2] )
