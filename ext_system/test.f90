program test
  use system_calls
  implicit none

  character(100) :: folder

  folder = "testfolder"

  call create_directory(trim(folder))

  open(101, file= trim(folder)//'/afile')
  write(101, *) 'test', 101
  close(101)

  !call remove_file(trim(folder)//'/afile')

  !call remove_directory(folder)


end program test
