program indi_star
!! Test program for individual stellar physics

    call do_main
    
    stop
end program

subroutine do_main
    use stellar_tracks_and_isochrones

    call test_routine

end subroutine do_main