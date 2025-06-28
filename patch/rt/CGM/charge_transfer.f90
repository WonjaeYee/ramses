!     Code by Jim Kingdon, in collaboration with G.J. Ferland
!     Modified by Harley Katz to a f90 module

MODULE charge_transfer
    use amr_parameters, only: dp
    implicit none

    private   ! default
    public init_ct_tables, HCTRecom, HCTIon, HEIIRecomb, HEIIIRecomb, C_HE_Recomb, C_HE_Ion &
        , N_HE_Recomb, N_HE_Ion, O_HE_Recomb, O_HE_Ion, MG_HE_Recomb      &
        , MG_HE_Ion, HEIIon

    real(dp), dimension(6,4,30), save:: CTRecomb = 0.d0
    real(dp), dimension(7,4,30), save:: CTIon = 0.d0

    CONTAINS
    SUBROUTINE init_ct_tables()
       implicit none

        CTRecomb = 0.d0
        !     Note: First parameter is in units of 1e-9!
        CTRecomb(:,1,2) = (/ 7.47e-6,2.06,9.93,-3.89,6e3,1e5 /)
        CTRecomb(:,2,2) = (/ 1.00e-5,0.,0.,0.,1e3,1e7 /)
        CTRecomb(:,1,3) = (/ 0.,0.,0.,0.,0.,0. /)
        CTRecomb(:,2,3) = (/ 1.26,0.96,3.02,-0.65,1e3,3e4 /)
        CTRecomb(:,3,3) = (/ 1.00e-5,0.,0.,0.,2e3,5e4 /)
        CTRecomb(:,1,4) = (/ 0.,0.,0.,0.,0.,0. /)
        CTRecomb(:,2,4) = (/ 1.00e-5,0.,0.,0.,2e3,5e4 /)
        CTRecomb(:,3,4) = (/ 1.00e-5,0.,0.,0.,2e3,5e4 /)
        CTRecomb(:,4,4) = (/ 5.17,0.82,-0.69,-1.12,2e3,5e4 /)
        CTRecomb(:,1,5) = (/ 0.,0.,0.,0.,0.,0. /)
        CTRecomb(:,2,5) = (/ 2.00e-2,0.,0.,0.,1e3,1e9 /)
        CTRecomb(:,3,5) = (/ 1.00e-5,0.,0.,0.,2e3,5e4 /)
        CTRecomb(:,4,5) = (/ 2.74,0.93,-0.61,-1.13,2e3,5e4 /)
        CTRecomb(:,1,6) = (/ 4.88e-7,3.25,-1.12,-0.21,5.5e3,1e5 /)
        CTRecomb(:,2,6) = (/ 1.67e-4,2.79,304.72,-4.07,5e3,5e4 /)
        CTRecomb(:,3,6) = (/ 3.25,0.21,0.19,-3.29,1e3,1e5 /)
        CTRecomb(:,4,6) = (/ 332.46,-0.11,-9.95e-1,-1.58e-3,1e1,1e5 /)
        CTRecomb(:,1,7) = (/ 1.01e-3,-0.29,-0.92,-8.38,1e2,5e4 /)
        CTRecomb(:,2,7) = (/ 3.05e-1,0.60,2.65,-0.93,1e3,1e5 /)
        CTRecomb(:,3,7) = (/ 4.54,0.57,-0.65,-0.89,1e1,1e5 /)
        CTRecomb(:,4,7) = (/ 2.95,0.55,-0.39,-1.07,1e3,1e6 /)
        CTRecomb(:,1,8) = (/ 1.04,3.15e-2,-0.61,-9.73,1e1,1e4 /)
        CTRecomb(:,2,8) = (/ 1.04,0.27,2.02,-5.92,1e2,1e5 /)
        CTRecomb(:,3,8) = (/ 3.98,0.26,0.56,-2.62,1e3,5e4 /)
        CTRecomb(:,4,8) = (/ 2.52e-1,0.63,2.08,-4.16,1e3,3e4 /)
        CTRecomb(:,1,9) = (/ 0.,0.,0.,0.,0.,0. /)
        CTRecomb(:,2,9) = (/ 1.00e-5,0.,0.,0.,2e3,5e4 /)
        CTRecomb(:,3,9) = (/ 9.86,0.29,-0.21,-1.15,2e3,5e4 /)
        CTRecomb(:,4,9) = (/ 7.15e-1,1.21,-0.70,-0.85,2e3,5e4 /)
        CTRecomb(:,1,10) = (/ 0.,0.,0.,0.,0.,0. /)
        CTRecomb(:,2,10) = (/ 1.00e-5,0.,0.,0.,5e3,5e4 /)
        CTRecomb(:,3,10) = (/ 14.73,4.52e-2,-0.84,-0.31,5e3,5e4 /)
        CTRecomb(:,4,10) = (/ 6.47,0.54,3.59,-5.22,1e3,3e4 /)
        CTRecomb(:,1,11) = (/ 0.,0.,0.,0.,0.,0. /)
        CTRecomb(:,2,11) = (/ 1.00e-5,0.,0.,0.,2e3,5e4 /)
        CTRecomb(:,3,11) = (/ 1.33,1.15,1.20,-0.32,2e3,5e4 /)
        CTRecomb(:,4,11) = (/ 1.01e-1,1.34,10.05,-6.41,2e3,5e4 /)
        CTRecomb(:,1,12) = (/ 0.,0.,0.,0.,0.,0. /)
        CTRecomb(:,2,12) = (/ 8.58e-5,2.49e-3,2.93e-2,-4.33,1e3,3e4 /)
        CTRecomb(:,3,12) = (/ 6.49,0.53,2.82,-7.63,1e3,3e4 /)
        CTRecomb(:,4,12) = (/ 6.36,0.55,3.86,-5.19,1e3,3e4 /)
        CTRecomb(:,1,13) = (/ 0.,0.,0.,0.,0.,0. /)
        CTRecomb(:,2,13) = (/ 1.00e-5,0.,0.,0.,1e3,3e4 /)
        CTRecomb(:,3,13) = (/ 7.11e-5,4.12,1.72e4,-22.24,1e3,3e4 /)
        CTRecomb(:,4,13) = (/ 7.52e-1,0.77,6.24,-5.67,1e3,3e4 /)
        CTRecomb(:,1,14) = (/ 0.,0.,0.,0.,0.,0. /)
        CTRecomb(:,2,14) = (/ 6.77,7.36e-2,-0.43,-0.11,5e2,1e5 /)
        CTRecomb(:,3,14) = (/ 4.90e-1,-8.74e-2,-0.36,-0.79,1e3,3e4 /)
        CTRecomb(:,4,14) = (/ 7.58,0.37,1.06,-4.09,1e3,5e4 /)
        CTRecomb(:,1,15) = (/ 0.,0.,0.,0.,0.,0. /)
        CTRecomb(:,2,15) = (/ 1.74e-4,3.84,36.06,-0.97,1e3,3e4 /)
        CTRecomb(:,3,15) = (/ 9.46e-2,-5.58e-2,0.77,-6.43,1e3,3e4 /)
        CTRecomb(:,4,15) = (/ 5.37,0.47,2.21,-8.52,1e3,3e4 /)
        CTRecomb(:,1,16) = (/ 3.82e-7,11.10,2.57e4,-8.22,1e3,1e4 /)
        CTRecomb(:,2,16) = (/ 1.00e-5,0.,0.,0.,1e3,3e4 /)
        CTRecomb(:,3,16) = (/ 2.29,4.02e-2,1.59,-6.06,1e3,3e4 /)
        CTRecomb(:,4,16) = (/ 6.44,0.13,2.69,-5.69,1e3,3e4 /)
        CTRecomb(:,1,17) = (/ 0.,0.,0.,0.,0.,0. /)
        CTRecomb(:,2,17) = (/ 1.00e-5,0.,0.,0.,1e3,3e4 /)
        CTRecomb(:,3,17) = (/ 1.88,0.32,1.77,-5.70,1e3,3e4 /)
        CTRecomb(:,4,17) = (/ 7.27,0.29,1.04,-10.14,1e3,3e4 /)
        CTRecomb(:,1,18) = (/ 0.,0.,0.,0.,0.,0. /)
        CTRecomb(:,2,18) = (/ 1.00e-5,0.,0.,0.,1e3,3e4 /)
        CTRecomb(:,3,18) = (/ 4.57,0.27,-0.18,-1.57,1e3,3e4 /)
        CTRecomb(:,4,18) = (/ 6.37,0.85,10.21,-6.22,1e3,3e4 /)
        CTRecomb(:,1,19) = (/ 0.,0.,0.,0.,0.,0. /)
        CTRecomb(:,2,19) = (/ 1.00e-5,0.,0.,0.,1e3,3e4 /)
        CTRecomb(:,3,19) = (/ 4.76,0.44,-0.56,-0.88,1e3,3e4 /)
        CTRecomb(:,4,19) = (/ 1.00e-5,0.,0.,0.,1e3,3e4 /)
        CTRecomb(:,1,20) = (/ 0.,0.,0.,0.,0.,0. /)
        CTRecomb(:,2,20) = (/ 0.,0.,0.,0.,1e1,1e9 /)
        CTRecomb(:,3,20) = (/ 3.17e-2,2.12,12.06,-0.40,1e3,3e4 /)
        CTRecomb(:,4,20) = (/ 2.68,0.69,-0.68,-4.47,1e3,3e4 /)
        CTRecomb(:,1,21) = (/ 0.,0.,0.,0.,0.,0. /)
        CTRecomb(:,2,21) = (/ 0.,0.,0.,0.,1e1,1e9 /)
        CTRecomb(:,3,21) = (/ 7.22e-3,2.34,411.50,-13.24,1e3,3e4 /)
        CTRecomb(:,4,21) = (/ 1.20e-1,1.48,4.00,-9.33,1e3,3e4 /)
        CTRecomb(:,1,22) = (/ 0.,0.,0.,0.,0.,0. /)
        CTRecomb(:,2,22) = (/ 0.,0.,0.,0.,1e1,1e9 /)
        CTRecomb(:,3,22) = (/ 6.34e-1,6.87e-3,0.18,-8.04,1e3,3e4 /)
        CTRecomb(:,4,22) = (/ 4.37e-3,1.25,40.02,-8.05,1e3,3e4 /)
        CTRecomb(:,1,23) = (/ 0.,0.,0.,0.,0.,0. /)
        CTRecomb(:,2,23) = (/ 1.00e-5,0.,0.,0.,1e3,3e4 /)
        CTRecomb(:,3,23) = (/ 5.12,-2.18e-2,-0.24,-0.83,1e3,3e4 /)
        CTRecomb(:,4,23) = (/ 1.96e-1,-8.53e-3,0.28,-6.46,1e3,3e4 /)
        CTRecomb(:,1,24) = (/ 0.,0.,0.,0.,0.,0. /)
        CTRecomb(:,2,24) = (/ 5.27e-1,0.61,-0.89,-3.56,1e3,3e4 /)
        CTRecomb(:,3,24) = (/ 10.90,0.24,0.26,-11.94,1e3,3e4 /)
        CTRecomb(:,4,24) = (/ 1.18,0.20,0.77,-7.09,1e3,3e4 /)
        CTRecomb(:,1,25) = (/ 0.,0.,0.,0.,0.,0. /)
        CTRecomb(:,2,25) = (/ 1.65e-1,6.80e-3,6.44e-2,-9.70,1e3,3e4 /)
        CTRecomb(:,3,25) = (/ 14.20,0.34,-0.41,-1.19,1e3,3e4 /)
        CTRecomb(:,4,25) = (/ 4.43e-1,0.91,10.76,-7.49,1e3,3e4 /)
        CTRecomb(:,1,26) = (/ 0.,0.,0.,0.,0.,0. /)
        CTRecomb(:,2,26) = (/ 1.26,7.72e-2,-0.41,-7.31,1e3,1e5 /)
        CTRecomb(:,3,26) = (/ 3.42,0.51,-2.06,-8.99,1e3,1e5 /)
        CTRecomb(:,4,26) = (/ 14.60,3.57e-2,-0.92,-0.37,1e3,3e4 /)
        CTRecomb(:,1,27) = (/ 0.,0.,0.,0.,0.,0. /)
        CTRecomb(:,2,27) = (/ 5.30,0.24,-0.91,-0.47,1e3,3e4 /)
        CTRecomb(:,3,27) = (/ 3.26,0.87,2.85,-9.23,1e3,3e4 /)
        CTRecomb(:,4,27) = (/ 1.03,0.58,-0.89,-0.66,1e3,3e4 /)
        CTRecomb(:,1,28) = (/ 0.,0.,0.,0.,0.,0. /)
        CTRecomb(:,2,28) = (/ 1.05,1.28,6.54,-1.81,1e3,1e5 /)
        CTRecomb(:,3,28) = (/ 9.73,0.35,0.90,-5.33,1e3,3e4 /)
        CTRecomb(:,4,28) = (/ 6.14,0.25,-0.91,-0.42,1e3,3e4 /)
        CTRecomb(:,1,29) = (/ 0.,0.,0.,0.,0.,0. /)
        CTRecomb(:,2,29) = (/ 1.47e-3,3.51,23.91,-0.93,1e3,3e4 /)
        CTRecomb(:,3,29) = (/ 9.26,0.37,0.40,-10.73,1e3,3e4 /)
        CTRecomb(:,4,29) = (/ 11.59,0.20,0.80,-6.62,1e3,3e4 /)
        CTRecomb(:,1,30) = (/ 0.,0.,0.,0.,0.,0. /) 
        CTRecomb(:,2,30) = (/ 1.00e-5,0.,0.,0.,1e3,3e4 /)
        CTRecomb(:,3,30) = (/ 6.96e-4,4.24,26.06,-1.24,1e3,3e4 /)
        CTRecomb(:,4,30) = (/ 1.33e-2,1.56,-0.92,-1.20,1e3,3e4 /)

        CTIon = 0.d0
        CTIon(:,1,3) = (/2.84e-3,1.99,375.54,-54.07,1e2,1e4,0.0 /)
        CTIon(:,2,3) = (/ 0.,0.,0.,0.,0.,0.,0. /)
        CTIon(:,3,3) = (/ 0.,0.,0.,0.,0.,0.,0. /)
        CTIon(:,1,4) = (/ 0.,0.,0.,0.,0.,0.,0. /)
        CTIon(:,2,4) = (/ 0.,0.,0.,0.,0.,0.,0. /)
        CTIon(:,3,4) = (/ 0.,0.,0.,0.,0.,0.,0. /)
        CTIon(:,1,5) = (/ 0.,0.,0.,0.,0.,0.,0. /)
        CTIon(:,2,5) = (/ 0.,0.,0.,0.,0.,0.,0. /)
        CTIon(:,3,5) = (/ 0.,0.,0.,0.,0.,0.,0. /)
        CTIon(:,1,6) = (/ 1.07e-6,3.15,176.43,-4.29,1e3,1e5,0.0 /)
        CTIon(:,2,6) = (/ 0.,0.,0.,0.,0.,0.,0. /)
        CTIon(:,3,6) = (/ 0.,0.,0.,0.,0.,0.,0. /)
        CTIon(:,1,7) = (/ 4.55e-3,-0.29,-0.92,-8.38,1e2,5e4,1.086 /)
        CTIon(:,2,7) = (/ 0.,0.,0.,0.,0.,0.,0. /)
        CTIon(:,3,7) = (/ 0.,0.,0.,0.,0.,0.,0. /)
        CTIon(:,1,8) = (/ 7.40e-2,0.47,24.37,-0.74,1e1,1e4,0.023 /)
        CTIon(:,2,8) = (/ 0.,0.,0.,0.,0.,0.,0. /)
        CTIon(:,3,8) = (/ 0.,0.,0.,0.,0.,0.,0. /)
        CTIon(:,1,9) = (/ 0.,0.,0.,0.,0.,0.,0. /)
        CTIon(:,2,9) = (/ 0.,0.,0.,0.,0.,0.,0. /)
        CTIon(:,3,9) = (/ 0.,0.,0.,0.,0.,0.,0. /)
        CTIon(:,1,10) = (/ 0.,0.,0.,0.,0.,0.,0. /)
        CTIon(:,2,10) = (/ 0.,0.,0.,0.,0.,0.,0. /)
        CTIon(:,3,10) = (/ 0.,0.,0.,0.,0.,0.,0. /)
        CTIon(:,1,11) = (/ 3.34e-6,9.31,2632.31,-3.04,1e3,2e4,0.0 /)
        CTIon(:,2,11) = (/ 0.,0.,0.,0.,0.,0.,0. /)
        CTIon(:,3,11) = (/ 0.,0.,0.,0.,0.,0.,0. /)
        CTIon(:,1,12) = (/ 9.76e-3,3.14,55.54,-1.12,5e3,3e4,0.0 /)
        CTIon(:,2,12) = (/ 7.60e-5,0.00,-1.97,-4.32,1e4,3e5,1.670 /)
        CTIon(:,3,12) = (/ 0.,0.,0.,0.,0.,0.,0. /)
        CTIon(:,1,13) = (/ 0.,0.,0.,0.,0.,0.,0. /)
        CTIon(:,2,13) = (/ 0.,0.,0.,0.,0.,0.,0. /)
        CTIon(:,3,13) = (/ 0.,0.,0.,0.,0.,0.,0. /)
        CTIon(:,1,14) = (/ 0.92,1.15,0.80,-0.24,1e3,2e5,0.0 /)
        CTIon(:,2,14) = (/ 2.26,7.36e-2,-0.43,-0.11,2e3,1e5,3.031 /)
        CTIon(:,3,14) = (/ 0.,0.,0.,0.,0.,0.,0. /)
        CTIon(:,1,15) = (/ 0.,0.,0.,0.,0.,0.,0. /)
        CTIon(:,2,15) = (/ 0.,0.,0.,0.,0.,0.,0. /)
        CTIon(:,3,15) = (/ 0.,0.,0.,0.,0.,0.,0. /)
        CTIon(:,1,16) = (/ 1.00e-5,0.00,0.00,0.00,1e3,1e4,0.0 /)
        CTIon(:,2,16) = (/ 0.,0.,0.,0.,0.,0.,0. /)
        CTIon(:,3,16) = (/ 0.,0.,0.,0.,0.,0.,0. /)
        CTIon(:,1,17) = (/ 0.,0.,0.,0.,0.,0.,0. /)
        CTIon(:,2,17) = (/ 0.,0.,0.,0.,0.,0.,0. /)
        CTIon(:,3,17) = (/ 0.,0.,0.,0.,0.,0.,0. /)
        CTIon(:,1,18) = (/ 0.,0.,0.,0.,0.,0.,0. /)
        CTIon(:,2,18) = (/ 0.,0.,0.,0.,0.,0.,0. /)
        CTIon(:,3,18) = (/ 0.,0.,0.,0.,0.,0.,0. /)
        CTIon(:,1,19) = (/ 0.,0.,0.,0.,0.,0.,0. /)
        CTIon(:,2,19) = (/ 0.,0.,0.,0.,0.,0.,0. /)
        CTIon(:,3,19) = (/ 0.,0.,0.,0.,0.,0.,0. /)
        CTIon(:,1,20) = (/ 0.,0.,0.,0.,0.,0.,0. /)
        CTIon(:,2,20) = (/ 0.,0.,0.,0.,0.,0.,0. /)
        CTIon(:,3,20) = (/ 0.,0.,0.,0.,0.,0.,0. /)
        CTIon(:,1,21) = (/ 0.,0.,0.,0.,0.,0.,0. /)
        CTIon(:,2,21) = (/ 0.,0.,0.,0.,0.,0.,0. /)
        CTIon(:,3,21) = (/ 0.,0.,0.,0.,0.,0.,0. /)
        CTIon(:,1,22) = (/ 0.,0.,0.,0.,0.,0.,0. /)
        CTIon(:,2,22) = (/ 0.,0.,0.,0.,0.,0.,0. /)
        CTIon(:,3,22) = (/ 0.,0.,0.,0.,0.,0.,0. /)
        CTIon(:,1,23) = (/ 0.,0.,0.,0.,0.,0.,0. /)
        CTIon(:,2,23) = (/ 0.,0.,0.,0.,0.,0.,0. /)
        CTIon(:,3,23) = (/ 0.,0.,0.,0.,0.,0.,0. /)
        CTIon(:,1,24) = (/ 0.,0.,0.,0.,0.,0.,0. /)
        CTIon(:,2,24) = (/ 4.39,0.61,-0.89,-3.56,1e3,3e4,3.349 /)
        CTIon(:,3,24) = (/ 0.,0.,0.,0.,0.,0.,0. /)
        CTIon(:,1,25) = (/ 0.,0.,0.,0.,0.,0.,0. /)
        CTIon(:,2,25) = (/ 2.83e-1,6.80e-3,6.44e-2,-9.70,1e3,3e4,2.368 /)
        CTIon(:,3,25) = (/ 0.,0.,0.,0.,0.,0.,0. /)
        CTIon(:,1,26) = (/ 0.,0.,0.,0.,0.,0.,0. /)
        CTIon(:,2,26) = (/ 2.10,7.72e-2,-0.41,-7.31,1e4,1e5,3.005 /)
        CTIon(:,3,26) = (/ 0.,0.,0.,0.,0.,0.,0. /)
        CTIon(:,1,27) = (/ 0.,0.,0.,0.,0.,0.,0. /)
        CTIon(:,2,27) = (/ 1.20e-2,3.49,24.41,-1.26,1e3,3e4,4.044 /)
        CTIon(:,3,27) = (/ 0.,0.,0.,0.,0.,0.,0. /)
        CTIon(:,1,28) = (/ 0.,0.,0.,0.,0.,0.,0. /)
        CTIon(:,2,28) = (/ 0.,0.,0.,0.,0.,0.,0. /)
        CTIon(:,3,28) = (/ 0.,0.,0.,0.,0.,0.,0. /)
        CTIon(:,1,29) = (/ 0.,0.,0.,0.,0.,0.,0. /)
        CTIon(:,2,29) = (/ 0.,0.,0.,0.,0.,0.,0. /)
        CTIon(:,3,29) = (/ 0.,0.,0.,0.,0.,0.,0. /)
        CTIon(:,1,30) = (/ 0.,0.,0.,0.,0.,0.,0. /)
        CTIon(:,2,30) = (/ 0.,0.,0.,0.,0.,0.,0. /)
        CTIon(:,3,30) = (/ 0.,0.,0.,0.,0.,0.,0. /)
    END SUBROUTINE init_ct_tables
    
    FUNCTION HCTRecom(ion,nelem,te)
        ! ion is stage of ionization, 2 for the ion going to the atom
        ! nelem is atomic number of element, 2 up to 30
        ! Example:  O+ + H => O + H+ is HCTRecom(2,8,4)
        ! Note that temperature is in linear scale
        implicit none
        integer, intent(in)::ion, nelem
        real(dp), intent(in)::te 
        real(dp)::HCTRecom

        !     local variables
        real(dp)::tused 
        integer::ipIon
        !real(dp), dimension(6,4,30):: CTRecomb
        real(dp)::a_op,b_op,c_op,d_op,e_op,f_op,logT

        HCTRecom = 0.d0

        ! No charge transfer above 10^5 K
        if (te.gt.1.d5) then
            return
        endif

        if (nelem.eq.8) then ! deal with oxygen separately            
            if (ion.eq.2) then 
                if (te.lt.10.d0) then
                    HCTRecom = 3.744e-10
                    return
                endif
                a_op = 2.3344302e-10
                b_op = 2.3651505e-10
                c_op = -1.3146803e-10
                d_op = 2.9979994e-11
                e_op = -2.8577012e-12
                f_op = 1.1963502e-13
                logT = log(te)
                HCTRecom = ((((f_op*logT + e_op)*logT + d_op)*logT + c_op)*logT + b_op)*logT + a_op
                return 
            else if (ion.eq.3) then
                if (te.le.1500.d0) then
                    HCTRecom = 0.5337d-9 * ((te/100.d0)**(-0.076d0))
                else
                    HCTRecom = 0.4344d-9 + 0.6340d-9 * (log10(te/1500.d0)**2.06d0)
    
                endif
                return
            endif
        endif

        ipIon = ion - 1

        ! use statistical charge transfer for ion > 4
        if( ipIon.gt.4 ) then
            HCTRecom = 1.92e-9 * ipIon
            return
        endif

        ! Make sure te is between temp. boundaries; set constant outside of range
        tused = max( te,CTRecomb(5,ipIon,nelem) )
        tused = min( tused , CTRecomb(6,ipIon,nelem) )
        tused = tused * 1e-4

        ! the interpolation equation
        HCTRecom = CTRecomb(1,ipIon,nelem) * 1e-9 * (tused**CTRecomb(2,ipIon,nelem)) * (1. + CTRecomb(3,ipIon,nelem) * exp(CTRecomb(4,ipIon,nelem)*tused) )

    END FUNCTION HCTRecom

    FUNCTION HCTIon(ion,nelem,te)
    ! ion is stage of ionization, 1 for atom
    ! nelem is atomic number of element, 2 up to 30
    ! Example:  O + H+ => O+ + H is HCTIon(1,8,1e4)
    ! Note that temperature is in linear scale
    implicit none
    integer, intent(in)::ion, nelem
    real(dp), intent(in)::te 
    real(dp)::HCTIon

    ! local variables
    real(dp)::tused 
    integer::ipIon
    !real(dp), dimension(7,4,30):: CTIon
    real(dp)::a,b,c
    real(dp)::a_o,b_o,c_o,d_o,e_o,f_o,g_o,logT

    HCTIon = 0.d0

    ! No charge transfer above 10^5 K
    if (te.gt.1.d5) then
        return
    endif


    if (nelem.eq.8.and.ion.eq.1) then ! deal with oxygen separately
        if (te.le.10) then
            HCTIon = 4.749e-20
        else if (te.gt.10.and.te.le.190.) then
            a = -21.134531
            b = -242.06831
            c = 84.761441
            HCTIon = exp(a + (b/te) + (c/(te*te)))
        else if (te.gt.190..and.te.le.200.) then
            HCTIon = 2.18733e-12*(te-190.0) + 1.85823e-10
        else 
            a_o = -7.6767404e-14
            b_o = -3.7282001e-13
            c_o = -1.488594e-12
            d_o = -3.6606214e-12 
            e_o = 2.0699463e-12
            f_o = -2.6139493e-13
            g_o = 1.1580844e-14
            logT = log(te)
            HCTIon = (((((g_o*logT + f_o)*logT + e_o)*logT + d_o)*logT + c_o)*logT + b_o)*logT + a_o
        endif

        return
    endif

    ipIon = ion
    if( ipIon.gt.3 ) then
      HCTIon = 0.d0
      return
    endif

    ! Make sure te is between temp. boundaries; set constant outside of range
    tused = max( te,CTIon(5,ipIon,nelem) )
    tused = min( tused , CTIon(6,ipIon,nelem) )
    tused = tused * 1e-4
    tused = max(tused,1.d-10) ! harley added to prevent zero temperature

    ! the interpolation equation
    HCTIon = CTIon(1,ipIon,nelem) * 1e-9 * (tused**CTIon(2,ipIon,nelem)) * (1. + CTIon(3,ipIon,nelem) * exp(CTIon(4,ipIon,nelem)*tused) ) * exp(-1.d0 * CTIon(7,ipIon,nelem)/tused)

    END FUNCTION HCTIon

    FUNCTION HEIIon(te)
        !Glover & Brand 2003, Kimura 1993 (from yoshida 2006)
        !He + H+ --> He+ + H
        real(dp), intent(in)::te 
        real(dp):: HEIIon

        HEIIon = 0.d0

        ! No charge transfer above 10^5 K
        if (te.gt.1.d5) then
            return
        endif
        ! This particular rate seems to get very small at low te
        if (te.lt.3.d3) then
            return
        endif

        if (te.le.1.d4) then 
            HEIIon = 1.26d-9*(te**(-0.75d0))*exp(-1.275d5/te)
        else
            HEIIon = 4.d-37*(te**(4.74d0))
        endif
    END FUNCTION HEIIon

    FUNCTION HEIIRecomb(te)
        !Stibbe & Tennyson 1999 (from yoshida 2006)
        !He+ + H --> He + H+
        real(dp), intent(in)::te 
        real(dp):: HEIIRecomb

        HEIIRecomb = 0.d0

        ! No charge transfer above 10^5 K
        if (te.gt.1.d5) then
            return
        endif

        HEIIRecomb = 1.20d-15*((te/3d2)**0.25d0)
    END FUNCTION HEIIRecomb

    FUNCTION HEIIIRecomb(te)
        !The process H + He++ -> He+ + H+
        real(dp), intent(in)::te 
        real(dp):: HEIIIRecomb

        HEIIIRecomb = 0.d0

        ! No charge transfer above 10^5 K
        !if (te.gt.1.d5) then
        !    return
        !endif

        !HEIIIRecomb = HCTRecom(3,2,te)
    END FUNCTION HEIIIRecomb

    ! CARBON AND HELIUM
    FUNCTION C_HE_Recomb(im,te)
        integer, intent(in)::im
        real(dp), intent(in)::te 
        real(dp):: C_HE_Recomb

        C_HE_Recomb = 0.d0

        if (im.eq.4) C_HE_Recomb = 4.6d-19 * te * te
        if (im.eq.5) C_HE_Recomb = 1d-14
    END FUNCTION C_HE_Recomb

    FUNCTION C_HE_Ion(im,te)
        integer, intent(in)::im
        real(dp), intent(in)::te 
        real(dp):: C_HE_Ion

        C_HE_Ion = 0.d0

        if (im.eq.1) C_HE_Ion = 6.3d-15 * ( (max(te,1.d9)/300.d0)**0.75d0 )
        ! harley edited out...not well behaved at either temeprature end
        !if (im.eq.2) C_HE_Ion = 5d-20 * te * te * exp(0.07d-4*te) * exp(6.29d0/(te/11604.525d0))
    END FUNCTION C_HE_Ion

    ! NITROGEN AND HELIUM
    FUNCTION N_HE_Recomb(im,te)
        integer, intent(in)::im
        real(dp), intent(in)::te 
        real(dp):: N_HE_Recomb

        N_HE_Recomb = 0.d0

        if (im.eq.3) N_HE_Recomb = 0.8d-10
        if (im.eq.4) N_HE_Recomb = 1.5d-10
        if (im.eq.5) N_HE_Recomb = 2d-9
    END FUNCTION N_HE_Recomb

    FUNCTION N_HE_Ion(im,te)
        integer, intent(in)::im
        real(dp), intent(in)::te 
        real(dp):: N_HE_Ion

        N_HE_Ion = 0.d0

        ! harley edited this out due to runaway at both temperature ends
        !if (im.eq.2) N_HE_Ion = 3.7d-20 * te * te * exp(0.063d-4*te) * exp(1.44d0/(te/11604.525d0))
    END FUNCTION N_HE_Ion

    ! OXYGEN AND HELIUM
    FUNCTION O_HE_Recomb(im,te)
        integer, intent(in)::im
        real(dp), intent(in)::te 
        real(dp):: O_HE_Recomb

        O_HE_Recomb = 0.d0

        if (im.eq.3) O_HE_Recomb = 3.2d-14 * max(te,1d5) / (te**0.05d0) ! harley inserted max to prevent runaway
        if (im.eq.4) O_HE_Recomb = 1.d-9
        if (im.eq.5) O_HE_Recomb = 6.d-10
    END FUNCTION O_HE_Recomb

    FUNCTION O_HE_Ion(im,te)
        integer, intent(in)::im
        real(dp), intent(in)::te 
        real(dp):: O_HE_Ion
        real(dp):: tloc

        O_HE_Ion = 0.d0

        if (im.eq.1) then 
            tloc = max(te,1.d5) ! Upper limit to prevent runaway
            O_HE_Ion = 4.991d-15 * ((tloc/1d4)**0.3794d0) * exp(tloc/1.121d6)
            O_HE_Ion = O_HE_Ion + 2.780d-15 * ((tloc/1d4)**-0.2163d0) * exp( -1.d0 * min(1e7,tloc)/(-8.158d5) )
        endif
    END FUNCTION O_HE_Ion

    ! MAGNESIUM AND HELIUM
    FUNCTION MG_HE_Recomb(im,te)
        integer, intent(in)::im
        real(dp), intent(in)::te 
        real(dp):: MG_HE_Recomb

        MG_HE_Recomb = 0.d0

        if (im.eq.4) MG_HE_Recomb = 7.5d-10
        if (im.eq.5) MG_HE_Recomb = 1.4e-10 * (te**0.3)
    END FUNCTION MG_HE_Recomb

    FUNCTION MG_HE_Ion(im,te)
        integer, intent(in)::im
        real(dp), intent(in)::te 
        real(dp):: MG_HE_Ion

        MG_HE_Ion = 0.d0
    END FUNCTION MG_HE_Ion

END MODULE charge_transfer
