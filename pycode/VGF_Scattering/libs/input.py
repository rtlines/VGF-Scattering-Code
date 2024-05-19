import logging

import numpy as np
import numpy.typing as npt

logger = logging.getLogger(__name__)

def read_input_file(workfile: str = "") -> tuple[int, str, float, float, float, float, float, float, float, float, float, npt.NDArray[np.double], npt.NDArray[np.double]]:
    # comments is the coment string from vgfin
    # number_of_dipoles is the number of dipoles to use
    # wavelength is the wavelength in relative uints (as long as all lenght units 
    #      are the same it doesn't matter, but use microns)
    # alpha, beta, and gamma are the euler rotation angles
    # psi is the polarization angle (I think)
    # semimajor_axis is the particle semimajor axis (use microns)
    # permitivity_real and permitivity_complex are the complex components of the particle permitivity
    # dipole_size is the dipole size
    # cell_locations is a number_of_dipoles by 3 array of cell locations
    # cell_weights is a number_of_dipoles array of cell weights
    #
    # Start reading in data
    #print("start read intput file")
    
    # Initialize all variables and give type hints
    comments: str = ""
    
    number_of_dipoles: int = 0
    
    wavelength: float = 0.

    alpha: float = 0.
    beta : float = 0.
    gamma: float = 0.
    
    psi: float = 0.

    semimajor_axis: float = 0.

    permitivity_real   : float = 0.
    permitivity_complex: float = 0.

    dipole_size: float = 0.

    cell_locations: npt.NDArray[np.double] = np.zeros((number_of_dipoles,3), dtype=np.double)
    cell_weights  : npt.NDArray[np.double] = np.zeros(number_of_dipoles, dtype=np.double)
    
    # Rudimentary error handling, so the function still runs if no workfile is found
    try:    
        # Open the file and begin reading in values
        logger.debug(f"Opening {workfile}...")
        with open(workfile, "r") as f:
            comments = f.readline().strip()            # read a whole line
            
            # Second line is how many dipoles the particle will be broken into
            number_of_dipoles = int(f.readline())
            # print("number_of_dipoles",number_of_dipoles)
            
            # Next line contains orientation information
            read_data: str = f.readline()
            split_string: list[str] = read_data.strip().split(" ")  # parse the line
            values: list[float] = [float(value) for value in split_string if value] # List comprehension to remove empty items from the list
            wavelength, alpha, beta, gamma, psi, semimajor_axis = values # Unpacks the values to their respective variables
     
            read_data: str = f.readline()
            split_string: list[str] = read_data.strip().split(" ")  # parse the line
            values: list[float] = [float(value) for value in split_string if value] # Convert non-empty strings to floats
            permitivity_real, permitivity_complex, dipole_size = values # Unpack the values into three variables
            
            # Replace placeholder arrays with real values
            logger.debug(f"Reading in {number_of_dipoles} cell locations and weights")
            cell_locations = np.zeros((number_of_dipoles,3))
            cell_weights = np.zeros(number_of_dipoles)
            for i in range(number_of_dipoles):             # loop to input all the theta and phi values
                read_data = f.readline()       # read in a whole line
                while "  " in read_data:
                    read_data = read_data.replace("  "," ")
                read_data.lstrip()           # remove leading space
                read_data.lstrip()           # remove leading space
                split_string = read_data.strip().split(" ") #parse the line splitting where
                                                    #there are spaces
                                                    #Each line starts with a space
                                                    #then a number, then a space
                                                    #then a vector and then a weight.
                                                    # So we will get
                                                    # several parts to our split
                                                    #a space, a string that is a 
                                                    #number, then a second string 
                                                    #that is a number, etc.
                values = [float(value) for value in split_string if value] # Convert non-empty strings to floats
                cell_locations[i][0] = values[0]
                cell_locations[i][1] = values[1]
                cell_locations[i][2] = values[2]
                cell_weights[i]      = values[3]
            logger.debug(f"Loaded cell locations and weights.")
        logger.debug(f"{workfile} closed.")
        
    except FileNotFoundError:
        logger.warning(f"File {workfile} not found! Using default values instead.")

    return number_of_dipoles, comments, wavelength, alpha, beta, gamma, psi, semimajor_axis, permitivity_real, permitivity_complex, dipole_size, cell_locations, cell_weights

def read_kv(workfile: str = "") -> tuple[int, npt.NDArray[np.double], npt.NDArray[np.double], float, float, int, int]:
    # Initialize all variables and provide type hints
    
    number_of_k_vectors: int = 0

    k_thetas: npt.NDArray[np.double] = np.zeros(number_of_k_vectors, dtype = np.double)
    k_phis  : npt.NDArray[np.double] = np.zeros(number_of_k_vectors, dtype = np.double)

    error     : float = 0.
    error_last: float = 0.

    m_count: int = 0 # What are m-count and k-count?
    k_count: int = 0

    # workfile ="Test_2_5_5.kv"
    try:
        logger.debug(f"Opening {workfile}...")
        with open(workfile, "r") as f:
            read_data: str = f.readline()
            number_of_k_vectors = int(read_data)

            k_thetas = np.zeros(number_of_k_vectors, dtype = np.double)
            k_phis   = np.zeros(number_of_k_vectors, dtype = np.double)

            # loop to input all the theta and phi values
            logger.debug(f"Loading in {number_of_k_vectors} k-vectors...")
            for i in range(number_of_k_vectors):             
                read_data = f.readline()
                # Parse string, removing newline character and converting non-empty strings to floats
                split_string = [float(value) for value in read_data.strip().split(" ") if value] 
                                                               
                theta=split_string[0]
                phi=split_string[1]

                k_thetas[i] = theta
                k_phis[i] = phi
            logger.debug(f"k-vectors loaded.")

            read_data = f.readline()       # The last line has dummy strings
                                           # for kth and kph but has at the end 
            split_string = read_data.strip().split(" ")
            values: list[str] = [value for value in split_string if value]
            error      = float(values[0])
            error_last = float(values[1])

            m_count = int(values[2])
            k_count = int(values[3])
        logger.debug(f"{workfile} closed.")

    except FileNotFoundError:
        logger.warning(f"File {workfile} not found! Using default values instead.")

    return number_of_k_vectors, k_thetas, k_phis, error, error_last, m_count, k_count

