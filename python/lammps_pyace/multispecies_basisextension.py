import logging
import numpy as np
import pickle
from importlib import resources
import re

from collections import defaultdict
from copy import deepcopy
from itertools import combinations, permutations, combinations_with_replacement, product
from typing import Dict, List, Union, Tuple

from lammps_pyace import ACEBBasisSet, BBasisConfiguration, BBasisFunctionSpecification, BBasisFunctionsSpecificationBlock

element_patt = re.compile("([A-Z][a-z]?)")

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

ALL = "ALL"
UNARY = "UNARY"
BINARY = "BINARY"
TERNARY = "TERNARY"
QUATERNARY = "QUATERNARY"
QUINARY = "QUINARY"
KEYWORDS = [ALL, UNARY, BINARY, TERNARY, QUATERNARY, QUINARY, 'number_of_functions_per_element']

NARY_MAP = {UNARY: 1, BINARY: 2, TERNARY: 3, QUATERNARY: 4, QUINARY: 5}
PERIODIC_ELEMENTS = chemical_symbols = [
    'H', 'He',
    'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
    'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar',
    'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
    'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',
    'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
    'In', 'Sn', 'Sb', 'Te', 'I', 'Xe',
    'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy',
    'Ho', 'Er', 'Tm', 'Yb', 'Lu',
    'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi',
    'Po', 'At', 'Rn',
    'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk',
    'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr',
    'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn', 'Nh', 'Fl', 'Mc',
    'Lv', 'Ts', 'Og']


def clean_bbasisconfig(initial_bbasisconfig):
    for block in initial_bbasisconfig.funcspecs_blocks:
        block.lmaxi = 0
        block.nradmaxi = 0
        block.nradbaseij = 0
        block.radcoefficients = []
        block.funcspecs = []


def reset_bbasisconfig(bconf):
    """set crad=delta_nk, func.coeffs=[0...]"""
    for block in bconf.funcspecs_blocks:
        block.set_func_coeffs(np.zeros_like(block.get_func_coeffs()))
        radcoefficients = np.array(block.radcoefficients)
        if len(radcoefficients.shape) == 3:
            # C_ nlk = delta_nk
            radcoefficients[:, :, :] = 0.0
            for nk in range(min(radcoefficients.shape[0], radcoefficients.shape[2])):
                radcoefficients[nk, :, nk] = 1.0
            block.set_radial_coeffs(radcoefficients.flatten())



def unify_to_minimized_indices(seq, shift=0):
    """
    Unify to minimized ordered sequence of indices
    """
    seq_map = {e: i + shift for i, e in enumerate(sorted(set(seq)))}
    return tuple([seq_map[e] for e in seq])


def unify_by_ordering(mus_comb, ns_comb):
    """
    Unify mus_comb and ns_comb to minimized-indices sequence, combine pairwise and  sort
    """
    return tuple(sorted(zip(unify_to_minimized_indices(mus_comb), unify_to_minimized_indices(ns_comb))))


def unify_absolute_by_ordering(mus_comb, ns_comb):
    """
    Unify mus_comb and ns_comb by combining pairwise and sort
    """
    return tuple(sorted(zip(mus_comb, ns_comb)))


def unify_mus_ns_comb(mus_comb, ns_comb):
    """
    Unify mus_comb, ns_comb by unifying to min-inds, combining, sorting and
    unifying the pair one more time to minimized-indices sequence
    """
    unif_comb = unify_by_ordering(mus_comb, ns_comb)
    return unify_to_minimized_indices(unif_comb)


def unify_absolute_mus_ns_comb(mus_comb, ns_comb):
    """
    Unify mus_comb, ns_comb by combining, sorting
    """
    unif_comb = unify_absolute_by_ordering(mus_comb, ns_comb)
    return unif_comb


def create_species_block_without_funcs(elements_vec: List[str], block_spec: Dict) -> BBasisFunctionsSpecificationBlock:
    """
    Create a BBasisFunctionsSpecificationBlock
    :param elements_vec: block's elements, i.e. (ele1, ele2, ...)
    :param block_spec: block specification dictionary
    :param crad_initializer: "delta" or "random"
    :return: BBasisFunctionsSpecificationBlock
    """
    block = BBasisFunctionsSpecificationBlock()

    block.block_name = " ".join(elements_vec)
    block.elements_vec = elements_vec

    # embedding
    if "fs_parameters" in block_spec:
        block.fs_parameters = block_spec["fs_parameters"]

    if "ndensity" in block_spec:
        block.ndensityi = block_spec["ndensity"]

    if "npot" in block_spec:
        block.npoti = block_spec["npot"]

    if "drho_core_cut" in block_spec:
        block.drho_cut = block_spec["drho_core_cut"]
    if "rho_core_cut" in block_spec:
        block.rho_cut = block_spec["rho_core_cut"]

    # bonds
    if "NameOfCutoffFunction" in block_spec:
        block.NameOfCutoffFunctionij = block_spec["NameOfCutoffFunction"]

    if "core-repulsion" in block_spec:
        block.core_rep_parameters = block_spec["core-repulsion"]

    if "rcut" in block_spec:
        block.rcutij = block_spec["rcut"]
    if "dcut" in block_spec:
        block.dcutij = block_spec["dcut"]
    if "radbase" in block_spec:
        block.radbase = block_spec["radbase"]
    if "radparameters" in block_spec:
        block.radparameters = block_spec["radparameters"]

    if "nradbase" in block_spec:
        block.nradbaseij = block_spec["nradbase"]
    if "nradmax" in block_spec:
        block.nradmaxi = block_spec["nradmax"]
    if "lmax" in block_spec:
        block.lmaxi = block_spec["lmax"]

    if "r_in" in block_spec:
        block.r_in = block_spec["r_in"]
        block.inner_cutoff_type = "distance"
    if "delta_in" in block_spec:
        block.delta_in = block_spec["delta_in"]
        block.inner_cutoff_type = "distance"
    if "inner_cutoff_type" in block_spec:
        block.inner_cutoff_type = block_spec["inner_cutoff_type"]
    else:
        block.inner_cutoff_type = "distance"

    # if crad_initializer == "delta":
    crad = np.zeros((block.nradmaxi, block.lmaxi + 1, block.nradbaseij))

    for n in range(0, min(block.nradmaxi, block.nradbaseij)):  # +1 is excluded, because ns=1...
        crad[n, :, n] = 1.0
    # elif crad_initializer == "random":
    #     crad = np.random.randn(block.nradmaxi, block.lmaxi + 1, block.nradbaseij)
    # else:
    #     raise ValueError("Unknown crad_initializer={}. Could be only 'delta' or 'random'".format(crad_initializer))

    block.radcoefficients = crad

    return block


def generate_species_keys(elements, r):

    keys = set()
    for el in elements:
        for comb in combinations_with_replacement(elements, r):
            keys.add((el,) + comb)
    return sorted(keys)


def generate_all_species_keys(elements):
    """
    Generate all ordered in [1:...] slice permutations of elements,
    that would be the species blocks names

    :param elements: list of elements (str)
    :return: list of all generated block names
    """
    keys = []
    nelements = len(elements)
    for r in range(0, nelements + 1):
        for key in generate_species_keys(elements, r):
            keys.append(key)
    return keys


def species_key_to_bonds(key):
    """
    Unify the tuple `key` to a list of bond pairs:

        A     -> [A,A]
        A:B   -> [A,B] [B,A],
        A:BC  -> [A,B] [B,A], [A,C], [C,A]
        A:BCD -> [A,B] [B,A], [A,C], [C,A], [A,D], [D,A]

    :param key: tuple of elements
    :return: list of bond pairs
    """
    if len(key) == 1:
        bonds = [(key[0], key[0])]
    else:
        k0 = key[0]
        rkeys = key[1:]
        bonds = []
        for rk in rkeys:
            bonds.append((k0, rk))
            bonds.append((rk, k0))
    return bonds


def create_multispecies_basis_config(potential_config: Dict,
                                     unif_mus_ns_to_lsLScomb_dict: Dict = None,
                                     func_coefs_initializer="zero",
                                     initial_basisconfig: BBasisConfiguration = None,
                                     overwrite_blocks_from_initial_bbasis=False
                                     ) -> BBasisConfiguration:

    overwrite_blocks_from_initial_bbasis = potential_config.get('overwrite_blocks_from_initial_bbasis',
                                                                overwrite_blocks_from_initial_bbasis)
                                                                
    if unif_mus_ns_to_lsLScomb_dict is None:
        with resources.files('lammps_pyace').joinpath('unif_mus_ns_to_lsLScomb_dict.pckl').open('rb') as f:
            unif_mus_ns_to_lsLScomb_dict = pickle.load(f)
            
    element_ndensity_dict = None
    if initial_basisconfig is not None:
        # extract embeddings from initial_basisconfig
        bas = ACEBBasisSet(initial_basisconfig)
        initial_embeddings = {}
        element_ndensity_dict = {}
        for el_ind, emb in bas.map_embedding_specifications.items():
            emb_dict = {}
            emb_dict["npot"] = emb.npoti
            emb_dict["ndensity"] = emb.ndensity
            emb_dict["fs_parameters"] = emb.FS_parameters
            emb_dict["rho_core_cut"] = emb.rho_core_cutoff
            emb_dict["drho_core_cut"] = emb.drho_core_cutoff
            initial_embeddings[bas.elements_name[el_ind]] = emb_dict
            element_ndensity_dict[bas.elements_name[el_ind]] = emb.ndensity

        if "embeddings" not in potential_config:
            potential_config["embeddings"] = {}

        embeddings = potential_config["embeddings"]
        for elm, emb_spec in initial_embeddings.items():
            if elm not in embeddings:
                embeddings[elm] = emb_spec
        potential_config["embeddings"] = embeddings

    blocks_list = create_multispecies_basisblocks_list(potential_config,
                                                       element_ndensity_dict=element_ndensity_dict,
                                                       func_coefs_initializer=func_coefs_initializer,
                                                       unif_mus_ns_to_lsLScomb_dict=unif_mus_ns_to_lsLScomb_dict,
                                                       )
    # compare with initial_basisconfig, if some blocks are missing in generated config - add them:
    if initial_basisconfig is not None:
        new_block_dict = {bl.block_name: bl for bl in blocks_list}
        if overwrite_blocks_from_initial_bbasis: # overwrite new blocks with old-ones
            for initial_block in initial_basisconfig.funcspecs_blocks:
                new_block_dict[initial_block.block_name] = initial_block
                log.info("Block {} is overwritten from initial potential".format(initial_block.block_name))
                if initial_block.block_name not in new_block_dict:
                    blocks_list.append(initial_block)
            blocks_list = [bl for bl in new_block_dict.values()]
        else: # only add missing blocks
            for initial_block in initial_basisconfig.funcspecs_blocks:
                if initial_block.block_name not in new_block_dict:
                    blocks_list.append(initial_block)
                    log.info("New block {} is added from initial potential".format(initial_block.block_name))

    new_basis_conf = BBasisConfiguration()
    new_basis_conf.deltaSplineBins = potential_config.get("deltaSplineBins", 0.001)
    new_basis_conf.funcspecs_blocks = blocks_list
    validate_bonds_nradmax_lmax_nradbase(new_basis_conf)
    new_basis_conf.validate(raise_exception=True)

    return new_basis_conf


def get_element_ndensity_dict(block_spec_dict):
    element_ndensity_dict = {}
    for el, spec_val in block_spec_dict.items():
        if len(el) == 1:
            element_ndensity_dict[el[0]] = spec_val['ndensity']
    return element_ndensity_dict


def generate_blocks_specifications_dict(potential_config: Dict) -> Dict:

    ### Embeddings
    ### possible keywords: ALL, UNARY  + el
    if "embeddings" in potential_config:
        embeddings_ext = generate_embeddings_ext(potential_config)
    else:
        embeddings_ext = {}
    ### Bonds
    ### possible keywords: ALL, UNARY, BINARY + (el,el)
    if "bonds" in potential_config:
        bonds_ext = generate_bonds_ext(potential_config)
    else:
        bonds_ext = {}
    ### Functions
    ### possible keywords: ALL, UNARY, BINARY, TERNARY, QUATERNARY, QUINARY + (el,el,...)
    if "functions" in potential_config:
        functions_ext = generate_functions_ext(potential_config)
    else:
        functions_ext = {}
    ### Update bonds specifications according to maximum observable nmax, lmax, nradbasemax in functions specifications
    bonds_ext = update_bonds_ext(bonds_ext, functions_ext)

    ### Combine together to have block_spec specs
    block_spec_dict = deepcopy(functions_ext)

    # update with embedding info
    for key, emb_ext_val in embeddings_ext.items():
        if key in block_spec_dict:
            block_spec_dict[key].update(emb_ext_val)
            
    # update with bond info
    for key, bonds_ext_val in bonds_ext.items():
        #if len(set(key)) == 1:
        #    key = (key[0],)
        if key in block_spec_dict:
            block_spec_dict[key].update(bonds_ext_val)
            
    return block_spec_dict


def generate_functions_ext(potential_config):
    elements = potential_config["elements"]
    elements = sorted(elements)

    functions = potential_config["functions"].copy()
    functions_ext = defaultdict(dict)

    if ALL in functions:
        max_rank = len(functions[ALL]['nradmax_by_orders'])
        for rank in range(0, max_rank+1):
            for species in generate_species_keys(elements, rank):
                if rank == 0:
                    functions_ext[species].update({})
                else:
                    for key in ['nradmax_by_orders', 'lmin_by_orders', 'lmax_by_orders']:
                        if key in functions[ALL]:
                            functions_ext[species][key] = functions[ALL][key][:rank]
            
    for nary_key, nary_val in NARY_MAP.items():
        if nary_key in functions:
            for key in generate_species_keys(elements, r=nary_val):
                functions_ext[key].update(functions[nary_key])
                
    for k in functions:
        if k not in KEYWORDS:
            if isinstance(k, str):  # single species string
                key = tuple(element_patt.findall(k))
            else:
                key = tuple(k)
            # TODO extend permutations
            functions_ext[key].update(functions[k])

    # drop all keys, that has no specifications
    #functions_ext = {k: v for k, v in functions_ext.items() if len(v) > 0}

    return functions_ext


def generate_bonds_ext(potential_config):
    elements = potential_config["elements"]
    elements = sorted(elements)

    bonds = potential_config["bonds"].copy()
    bonds_ext = {pair: {} for pair in product(elements, repeat=2)}
    if ALL in bonds:
        for pair in bonds_ext:
            bonds_ext[pair].update(bonds[ALL])
    if UNARY in bonds:
        for el in elements:
            bonds_ext[(el, el)].update(bonds[UNARY])
    if BINARY in bonds:
        for pair in permutations(elements, 2):
            bonds_ext[pair].update(bonds[BINARY])
    for pair in bonds:
        if pair not in KEYWORDS:  # assume that pair is valid (el1, el2)
            if isinstance(pair, str):
                # dpair= (pair, pair)
                # use regex to
                dpair = tuple(element_patt.findall(pair))
                if len(dpair) == 1:
                    dpair = (dpair[0], dpair[0])
            else:
                dpair = pair
            bonds_ext[dpair].update(bonds[pair])
            r_pair = tuple(reversed(dpair))
            bonds_ext[r_pair].update(bonds[pair])
    # drop all keys, that has no specifications
    bonds_ext = {k: v for k, v in bonds_ext.items() if len(v) > 0}
    return bonds_ext


def generate_embeddings_ext(potential_config):
    elements = potential_config["elements"]
    elements = sorted(elements)

    embeddings = potential_config["embeddings"].copy()
    embeddings_ext = {(el,): {} for el in elements}
    # ALL and UNARY behave identically
    if ALL in embeddings:
        for el in elements:
            embeddings_ext[(el,)].update(embeddings[ALL])
    if UNARY in embeddings:
        for el in elements:
            embeddings_ext[(el,)].update(embeddings[UNARY])
    for el, val in embeddings.items():
        if el in elements:
            embeddings_ext[(el,)].update(val)
        elif el not in [ALL, UNARY]:
            raise ValueError(f"{el} is not in specified elements: {elements}")

    # drop all keys, that has no specifications
    embeddings_ext = {k: v for k, v in embeddings_ext.items() if len(v) > 0}
    return embeddings_ext


def update_bonds_ext(bonds_ext, functions_ext):
    # if bonds_ext is empty - return it as it is
    if not bonds_ext:
        return bonds_ext

    bonds_ext_updated = deepcopy(bonds_ext)
    # run through functions specifications and update/validate bond's nradbase, nradmax, lmax
    for key, funcs_spec in functions_ext.items():
    
        if len(key) < 2: continue
    
        nradbasemax = max(funcs_spec['nradmax_by_orders'][:1])
        if len(funcs_spec['nradmax_by_orders'][1:]) > 0:
            nradmax = max(funcs_spec['nradmax_by_orders'][1:])
        else:
            nradmax = 0
        lmax = max(funcs_spec['lmax_by_orders'])
        
        if len(key) > 2:
            funcs_spec['nradmax'] = max(nradmax, nradbasemax)
            funcs_spec['nradbasemax'] = nradbasemax
            funcs_spec['lmax'] = lmax
            funcs_spec['nradbase'] = nradbasemax
            funcs_spec['rcut'] = np.inf
            funcs_spec['dcut'] = 0.0
            funcs_spec['rcut_in'] = 0.0
            funcs_spec['dcut_in'] = 0.0

        for bkey in species_key_to_bonds(key):
        
            bond = bonds_ext[bkey]
                    
            if len(key)>2:
                bonds_bkey = bonds_ext_updated[bkey]
                funcs_spec['radparameters'] = bonds_bkey['radparameters']
                funcs_spec['core-repulsion'] = bonds_bkey['core-repulsion']
                funcs_spec['rcut'] = min(funcs_spec['rcut'], bonds_bkey['rcut'])
                funcs_spec['dcut'] = max(funcs_spec['dcut'], bonds_bkey['dcut'])
                funcs_spec['rcut_in'] = max(funcs_spec['rcut_in'], bonds_bkey['rcut_in'])
                funcs_spec['dcut_in'] = max(funcs_spec['dcut_in'], bonds_bkey['dcut_in'])

            if 'nradbase' not in bond:
                if bonds_ext_updated[bkey].get('nradbase', 0) < nradbasemax:
                    bonds_ext_updated[bkey]['nradbase'] = nradbasemax
            else:
                if bond['nradbase'] < nradbasemax:
                    raise ValueError(f"nradbase={bond['nradbase']} for bond {bkey} " + \
                                     f"is less than nradbasemax={nradbasemax} from {key}")

            if 'nradmax' not in bond:
                if bonds_ext_updated[bkey].get('nradmax', 0) < nradmax:
                    bonds_ext_updated[bkey]['nradmax'] = nradmax
            else:
                if bond['nradmax'] < nradmax:
                    raise ValueError(f"nradmax={bond['nradmax']} for bond {bkey} " + \
                                     f"is less than nradmax={nradmax} from {key}")

            if 'lmax' not in bond:
                if bonds_ext_updated[bkey].get('lmax', 0) < lmax:
                    bonds_ext_updated[bkey]['lmax'] = lmax
            else:
                if bond['lmax'] < lmax:
                    raise ValueError(f"lmax={bond['lmax']} for bond {bkey} " + \
                                     f"is less than nradmax={lmax} from {key}")
                                     
    return bonds_ext_updated


def create_multispecies_basisblocks_list(potential_config: Dict,
                                         element_ndensity_dict: Dict = None,
                                         func_coefs_initializer="zero",
                                         unif_mus_ns_to_lsLScomb_dict=None,
                                         verbose=False) -> List[BBasisFunctionsSpecificationBlock]:
                                         
    blocks_specifications_dict = generate_blocks_specifications_dict(potential_config)
    
    #print(f"*** blocks_specifications_dict")
    #for k,v in blocks_specifications_dict.items():
    #    print(f"*** {k}: {v}")
    #print(f"")

    element_ndensity_dict =  element_ndensity_dict or {}
    constr_element_ndensity_dict = get_element_ndensity_dict(blocks_specifications_dict)
    for k, v in constr_element_ndensity_dict.items():
        if k not in element_ndensity_dict:
            element_ndensity_dict[k] = v
    if not element_ndensity_dict:
        raise ValueError("`element_ndensity_dict` neither provided nor constructed")

    blocks_list = []
    for elements_vec, block_spec_dict in blocks_specifications_dict.items():
        if verbose:
            print("Block elements:", elements_vec)
            
        ndensity = element_ndensity_dict[elements_vec[0]]
        spec_block = create_species_block(elements_vec, block_spec_dict, ndensity,
                                          func_coefs_initializer, unif_mus_ns_to_lsLScomb_dict)

        if verbose:
            print(len(spec_block.funcspecs), " functions added")
        blocks_list.append(spec_block)
    return blocks_list


def create_species_block(elements_vec: List, block_spec_dict: Dict,
                         ndensity: int,
                         func_coefs_initializer="zero",
                         unif_mus_ns_to_lsLScomb_dict=None) -> BBasisFunctionsSpecificationBlock:
    central_atom = elements_vec[0]

    elms = tuple(sorted(set(elements_vec)))
    nary = len(elms)
    spec_block = create_species_block_without_funcs(elements_vec, block_spec_dict)
    current_block_func_spec_list = []
    if "nradmax_by_orders" in block_spec_dict and "lmax_by_orders" in block_spec_dict:
        rank = len(elements_vec)-1
        unif_abs_combs_set = set()
        
        nmax = block_spec_dict["nradmax_by_orders"][rank-1]
        lmax = block_spec_dict["lmax_by_orders"][rank-1]
        
        if "lmin_by_orders" in block_spec_dict:
            lmin = block_spec_dict["lmin_by_orders"][rank-1]
        else:
            lmin = 0
            
        ns_range = range(1, nmax + 1)
        mus_comb = elements_vec[1:]
        mus_comb_ext = tuple([central_atom] + list(mus_comb))  # central atom + ordered tail
       
        for ns_comb in product(ns_range, repeat=rank):  # exhaustive list
            unif_abs_comb = unify_absolute_mus_ns_comb(mus_comb, ns_comb)
            if unif_abs_comb in unif_abs_combs_set:
                continue
            unif_abs_combs_set.add(unif_abs_comb)
            unif_comb = unify_mus_ns_comb(mus_comb, ns_comb)
            if unif_comb not in unif_mus_ns_to_lsLScomb_dict:
                raise ValueError(
                    "Specified potential shape is too big " + \
                    "and goes beyond the precomputed BBasisFunc white-list" + \
                    "for unified combination {}".format(unif_comb))

            mus_ns_white_list = unif_mus_ns_to_lsLScomb_dict[unif_comb]  # only ls, LS are important
                
            for (pre_ls, pre_LS) in mus_ns_white_list:
                if lmin <= min(pre_ls) and max(pre_ls) <= lmax:
                    if "coefs_init" in block_spec_dict:
                        func_coefs_initializer = block_spec_dict["coefs_init"]

                    if func_coefs_initializer == "zero":
                        coefs = [0] * ndensity
                    elif func_coefs_initializer == "random":
                        coefs = np.random.randn(ndensity) * 1e-4
                    else:
                        raise ValueError(
                            "Unknown func_coefs_initializer={}. Could be only 'zero' or 'random'".format(func_coefs_initializer))
                    
                    new_spec = BBasisFunctionSpecification(elements=mus_comb_ext,
                                                               ns=ns_comb,
                                                               ls=pre_ls,
                                                               LS=pre_LS,
                                                               coeffs=coefs
                                                            )

                    current_block_func_spec_list.append(new_spec)
                        
        spec_block.funcspecs = current_block_func_spec_list
    return spec_block























































def is_mult_basisfunc_equivalent(func1: BBasisFunctionSpecification, func2: BBasisFunctionSpecification) -> bool:
    return (func1.elements == func2.elements) and \
           (func1.ns == func2.ns) and \
           (func1.ls == func2.ls) and \
           (func1.LS == func2.LS)


class BlockBasisFunctionsList:

    def __init__(self, block: BBasisFunctionsSpecificationBlock):
        self.funcs = block.funcspecs

    def find_existing(self, func: BBasisFunctionSpecification) -> bool:
        for other_func in self.funcs:
            if is_mult_basisfunc_equivalent(func, other_func):
                return other_func
        return None


def extend_basis_block(init_block: BBasisFunctionsSpecificationBlock,
                       final_block: BBasisFunctionsSpecificationBlock,
                       num_funcs=None,
                       ladder_type="body_order") -> Tuple[BBasisFunctionsSpecificationBlock, bool]:
    # check that block name is identical
    assert init_block.block_name == final_block.block_name, ValueError("Could not extend block '' to new block ''". \
                                                                       format(init_block.block_name,
                                                                              final_block.block_name
                                                                              ))

    nelements = init_block.number_of_species

    init_block_list = BlockBasisFunctionsList(init_block)
    final_basis_funcs = sort_funcspecs_list(final_block.funcspecs, ladder_type)

    # existing_funcs_list = []
    existing_funcs_list = init_block.funcspecs

    new_funcs_list = []
    for new_func in final_basis_funcs:
        existing_func = init_block_list.find_existing(new_func)
        if existing_func is None:
            #             print("Func ", new_func, " added")
            new_funcs_list.append(new_func)

    if num_funcs is not None and len(new_funcs_list) > num_funcs:
        new_funcs_list = new_funcs_list[:num_funcs]
    else:
        new_funcs_list = new_funcs_list  # use all new funcs

    extended_block = init_block.copy()

    # if no new functions to add, return init_block
    if len(new_funcs_list) == 0:
        return extended_block, False

    extended_func_list = sort_funcspecs_list(existing_funcs_list + new_funcs_list, "body_order")

    # Update crad only for nelements<=2
    if nelements <= 2:
        validate_radial_shape_from_funcs(extended_block, extended_func_list)
        initialize_block_crad(extended_block)

        extended_radcoeffs = np.array(extended_block.radcoefficients)
        init_radcoeffs = np.array(init_block.radcoefficients)
        merge_crad_matrix(extended_radcoeffs, init_radcoeffs)

        extended_block.radcoefficients = extended_radcoeffs

    extended_block.funcspecs = extended_func_list

    # core-repulsion translating from final_basis

    if nelements <= 2:
        extended_block.core_rep_parameters = final_block.core_rep_parameters

    if nelements == 1:
        extended_block.rho_cut = final_block.rho_cut
        extended_block.drho_cut = final_block.drho_cut

    return extended_block, True


def merge_crad_matrix(extended_radcoeffs, init_radcoeffs):
    init_radcoeffs = np.array(init_radcoeffs)
    if len(init_radcoeffs.shape) == 3:
        common_shape = [min(s1, s2) for s1, s2 in zip(np.shape(init_radcoeffs), np.shape(extended_radcoeffs))]
        if len(common_shape) == 3:
            extended_radcoeffs[:common_shape[0], :common_shape[1], :common_shape[2]] = \
                init_radcoeffs[:common_shape[0], :common_shape[1], :common_shape[2]]


def initialize_block_crad(extended_block, crad_init="delta"):
    init_radcoeffs = np.array(extended_block.radcoefficients)
    if len(init_radcoeffs.shape) == 3:
        new_nradmax = max(extended_block.nradmaxi, init_radcoeffs.shape[0])
        new_lmax = max(extended_block.lmaxi, init_radcoeffs.shape[1] - 1)
        new_nradbase = max(extended_block.nradbaseij, init_radcoeffs.shape[2])
    else:
        new_nradbase = extended_block.nradbaseij
        new_lmax = extended_block.lmaxi
        new_nradmax = extended_block.nradmaxi

    if crad_init == "delta":
        extended_radcoeffs = np.zeros((new_nradmax, new_lmax + 1, new_nradbase))
        for n in range(min(new_nradmax, new_nradbase)):
            extended_radcoeffs[n, :, n] = 1.
    elif crad_init == "zero":
        extended_radcoeffs = np.zeros((new_nradmax, new_lmax + 1, new_nradbase))
    elif crad_init == "random":
        extended_radcoeffs = np.random.randn(*(new_nradmax, new_lmax + 1, new_nradbase))
    else:
        raise ValueError("Unknown value for crad_init ({}). Use delta, zero or random".format(crad_init))

    merge_crad_matrix(extended_radcoeffs, init_radcoeffs)

    extended_block.nradmaxi = new_nradmax
    extended_block.lmaxi = new_lmax
    extended_block.nradbaseij = new_nradbase

    extended_block.radcoefficients = extended_radcoeffs


def validate_radial_shape_from_funcs(extended_block, func_list=None):
    new_nradmax = 0
    new_nradbase = 0
    new_lmax = 0
    if func_list is None:
        func_list = extended_block.funcspecs
    for func in func_list:
        rank = len(func.ns)
        if rank == 1:
            new_nradbase = max(max(func.ns), new_nradbase)
        else:
            new_nradmax = max(max(func.ns), new_nradmax)
        new_lmax = max(max(func.ls), new_lmax)
    extended_block.nradbaseij = new_nradbase
    extended_block.lmaxi = new_lmax
    extended_block.nradmaxi = new_nradmax


def validate_bonds_nradmax_lmax_nradbase(ext_basis: BBasisConfiguration):
    """
    Check the nradbase, lmax and nradmax over all bonds in all species blocks
    """
    ext_blocks_dict = {block.block_name: block for block in ext_basis.funcspecs_blocks}

    max_nlk_dict = defaultdict(lambda: defaultdict(int))

    for block_name, block in ext_blocks_dict.items():
    
        # print(f"*** block_name {block_name} block {block}")

        for f in block.funcspecs:
            rank = len(f.ns)
            mu0 = f.elements[0]
            mus = f.elements[1:]
            ns = f.ns
            ls = f.ls

            for mu, n, l in zip(mus, ns, ls):
                bond = (mu0, mu)

                if rank == 1:
                    max_nlk_dict[bond]["nradbase"] = max(max_nlk_dict[bond]["nradbase"], n)
                else:
                    max_nlk_dict[bond]["nradmax"] = max(max_nlk_dict[bond]["nradmax"], n)

                max_nlk_dict[bond]["lmax"] = max(max_nlk_dict[bond]["lmax"], l)

    # loop over max_nlk_dict and symmetrize pair bonds
    for bond_pair, dct in max_nlk_dict.items():
        if len(bond_pair) == 2:
            sym_bond_pair = (bond_pair[1], bond_pair[0])
            sym_dct = max_nlk_dict[sym_bond_pair]
            max_nradbase = max(dct["nradbase"], sym_dct["nradbase"])
            max_lmax = max(dct["lmax"], sym_dct["lmax"])
            max_nradmax = max(dct["nradmax"], sym_dct["nradmax"])

            max_nlk_dict[bond_pair]["nradbase"] = max_nlk_dict[sym_bond_pair]["nradbase"] = max_nradbase
            max_nlk_dict[bond_pair]["lmax"] = max_nlk_dict[sym_bond_pair]["lmax"] = max_lmax
            max_nlk_dict[bond_pair]["nradmax"] = max_nlk_dict[sym_bond_pair]["nradmax"] = max_nradmax

    bonds_dict = {}
    for k, v in max_nlk_dict.items():
        #if k[0] == k[1]:
        #    bonds_dict[k[0]] = v
        #else:
        bonds_dict[" ".join(k)] = v
            
    for block in ext_basis.funcspecs_blocks:
        k = block.block_name
        
        if len(k.split()) != 2:
            continue

        # skip more ternary and higher blocks, because bond specification are defined only in unary/binary blocks
        if len(k.split()) > 2:
            #print(f"*** k {k} block {block}")
            #block.radcoefficients = []
            #block.nradbaseij = 0
            #block.lmaxi = 0
            #block.nradmaxi = 0
            continue

        if block.nradbaseij < bonds_dict[k]["nradbase"]:
            block.nradbaseij = bonds_dict[k]["nradbase"]

        if block.nradmaxi < bonds_dict[k]["nradmax"]:
            block.nradmaxi = bonds_dict[k]["nradmax"]

        if block.lmaxi < bonds_dict[k]["lmax"]:
            block.lmaxi = bonds_dict[k]["lmax"]

        initialize_block_crad(block)












