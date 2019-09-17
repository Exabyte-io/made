
//package burai.com.consts; }

export class Constants {

    /*
     * Physical constants, SI (NIST CODATA 2006), Web Version 5.1
     * <http://physics.nist.gov/constants>
     */
    static get H_PLANCK_SI() { return 6.62606896E-34; } // J s
    static get K_BOLTZMANN_SI() { return 1.3806504E-23; } // J K^-1
    static get ELECTRON_SI() { return 1.602176487E-19; } // C
    static get ELECTRONVOLT_SI() { return 1.602176487E-19; } // J
    static get ELECTRONMASS_SI() { return 9.10938215E-31; } // Kg
    static get HARTREE_SI() { return 4.35974394E-18; } // J
    static get RYDBERG_SI() { return Constants.HARTREE_SI / 2.0; } /* J */
    static get BOHR_RADIUS_SI() { return 0.52917720859E-10; } // m
    static get AMU_SI() { return 1.660538782E-27; } // Kg
    static get C_SI() { return 2.99792458E+8; } // m sec^-1
    static get MUNOUGHT_SI() { return 4.0 * Math.PI * 1.0E-7; } // N A^-2
    static get EPSNOUGHT_SI() { return 1.0 / (Constants.MUNOUGHT_SI * C_SI * C_SI); } /* F m^-1 */

    /*
     *
     * Physical constants, atomic units:
     * AU for "Hartree" atomic units (e() { return m() { return hbar() { return 1)
     * RY for "Rydberg" atomic units (e^2=2, m=1/2, hbar=1)
     */
    static get K_BOLTZMANN_AU() { return Constants.K_BOLTZMANN_SI / Constants.HARTREE_SI; }
    static get K_BOLTZMANN_RY() { return Constants.K_BOLTZMANN_SI / Constants.RYDBERG_SI; }

    /*
     * Unit conversion factors: energy and masses
     */
    static get AUTOEV() { return Constants.HARTREE_SI / Constants.ELECTRONVOLT_SI; }
    static get RYTOEV() { return Constants.AUTOEV / 2.0; }
    static get AMU_AU() { return Constants.AMU_SI / Constants.ELECTRONMASS_SI; }
    static get AMU_RY() { return Constants.AMU_AU / 2.0; }

    /*
     * Unit conversion factors: atomic unit of time, in s and ps
     */
    static get AU_SEC() { return Constants.H_PLANCK_SI / (2.0 * Math.PI) / Constants.HARTREE_SI; }
    static get AU_PS() { return Constants.AU_SEC * 1.0E+12; }

    /*
     * Unit conversion factors: pressure (1 Pa() { return 1 J/m^3, 1GPa() { return 10 Kbar )
     */
    static get AU_GPA() { return Constants.HARTREE_SI / Constants.BOHR_RADIUS_SI / Constants.BOHR_RADIUS_SI / Constants.BOHR_RADIUS_SI / 1.0E+9; }
    static get RY_KBAR() { return 10.0 * Constants.AU_GPA / 2.0; }

    /*
     * Unit conversion factors: 1 debye() { return 10^-18 esu*cm
     *                                 () { return 3.3356409519*10^-30 C*m
     *                                 () { return 0.208194346 e*A
     * ( 1 esu() { return (0.1/c) Am, c=299792458 m/s)
     */
    static get DEBYE_SI() { return 3.3356409519 * 1.0E-30; } // C*m
    static get AU_DEBYE() { return Constants.ELECTRON_SI * Constants.BOHR_RADIUS_SI / DEBYE_SI; }

    static get eV_to_kelvin() { return Constants.ELECTRONVOLT_SI / Constants.K_BOLTZMANN_SI; }
    static get ry_to_kelvin() { return Constants.RYDBERG_SI / Constants.K_BOLTZMANN_SI; }

    /*
     * Unit conversion factors: Energy to wavelength
     */
    static get EVTONM() { return 1.0E+9 * Constants.H_PLANCK_SI * Constants.C_SI / Constants.ELECTRONVOLT_SI; }
    static get RYTONM() { return 1.0E+9 * Constants.H_PLANCK_SI * Constants.C_SI / Constants.RYDBERG_SI; }

    /*
     * Speed of light in atomic units
     */
    static get C_AU() { return Constants.C_SI / Constants.BOHR_RADIUS_SI * Constants.AU_SEC; }

    /*
     * Temperature
     */
    static get CENTIGRADE_ZERO() { return 273.15; }

    /*
     * COMPATIBIILITY
     */
    static get BOHR_RADIUS_CM() { return Constants.BOHR_RADIUS_SI * 100.0; }
    static get BOHR_RADIUS_ANGS() { return Constants.BOHR_RADIUS_CM * 1.0E+8; }
    static get ANGSTROM_AU() { return 1.0 / Constants.BOHR_RADIUS_ANGS; }
    static get DIP_DEBYE() { return Constants.AU_DEBYE; }
    static get AU_TERAHERTZ() { return Constants.AU_PS; }
    static get AU_TO_OHMCMM1() { return 46000.0; } // (ohm cm)^-1
    static get RY_TO_THZ() { return 1.0 / Constants.AU_TERAHERTZ / (4.0 * Math.PI); }
    static get RY_TO_GHZ() { return Constants.RY_TO_THZ * 1000.0; }
    static get RY_TO_CMM1() { return 1.E+10 * Constants.RY_TO_THZ / Constants.C_SI; }
}
