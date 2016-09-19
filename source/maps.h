#define MAP_LINELENGTH 160
#define MAX_MAPS 10

// these definitions are for various ionization modes
#define IONMODE_ML93_FIXTE 0                  // Lucy Mazzali using existing t_e (no HC balance)
#define IONMODE_LTE 1                          // LTE using t_r
#define IONMODE_FIXED 2                       // Hardwired concentrations
#define IONMODE_ML93 3                        // Lucy Mazzali
#define IONMODE_LTE_SIM 4                     // LTE with SIM correction 
#define IONMODE_PAIRWISE_ML93 6               // pairwise version of Lucy Mazzali 
#define IONMODE_PAIRWISE_SPECTRALMODEL 7      // pariwise modeled J_nu approach
#define IONMODE_MATRIX_BB 8	              // matrix solver BB model
#define IONMODE_MATRIX_SPECTRALMODEL 9        // matrix solver spectral model

// and the corresponding modes in nebular_concentrations
#define NEBULARMODE_TR 0                       // LTE using t_r
#define NEBULARMODE_TE 1                       // LTE using t_e
#define NEBULARMODE_ML93 2                     // ML93 using correction
#define NEBULARMODE_NLTE_SIM 3                 // Non_LTE with SS modification (Probably could be removed)
#define NEBULARMODE_LTE_GROUND 4               // A test mode which foces all levels to the GS (Probably could be removed)
#define NEBULARMODE_PAIRWISE_ML93 6            // pairwise ML93
#define NEBULARMODE_PAIRWISE_SPECTRALMODEL 7   // pairwise spectral models
#define NEBULARMODE_MATRIX_BB 8	               // matrix solver BB model
#define NEBULARMODE_MATRIX_SPECTRALMODEL 9     // matrix solver spectral model

// arrays to use for ionization modes
struct mode_maps
{
  char ion_modes[MAX_MAPS][MAP_LINELENGTH];
}
maps;






