    program main
    implicit none 
    integer , save :: flag = 0
    double precision :: ran 
    if(flag==0) then 
      call random_seed()
      flag = 1 
    endif 
    call random_number(ran)     ! built in fortran 90 random number function 
     !!! returen number between 0 and 1
    write(*,*) ran
    
    end program main