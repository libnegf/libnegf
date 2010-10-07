program create_H

  integer :: i,j


  open(101,file='H_real.dat')
  open(102,file='H_imm.dat')  

  write(101,*) "% Size = 10 10" 
  write(101,*) "% Nonzeros = 18" 
  write(101,*) "zzz = zeros(18,3);"
  write(101,*) "zzz = ["

  write(102,*) "% Size = 10 10" 
  write(102,*) "% Nonzeros = 18" 
  write(102,*) "zzz = zeros(18,3);"
  write(102,*) "zzz = ["

  do i=1,9
     write(101,*) i,i+1,1.d0
     write(101,*) i+1,i,1.d0
     write(102,*) i,i+1,0.d0
     write(102,*) i+1,i,0.d0
  end do

  write(101,*) "]"
  write(102,*) "]"

end program create_H
