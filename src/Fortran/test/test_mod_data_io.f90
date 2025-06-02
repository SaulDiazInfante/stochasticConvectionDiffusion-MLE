program test_mod_data_io
    use mod_data_io
    implicit none
    
    call test_write_data_creates_file()

contains
    
    subroutine test_write_data_creates_file()
        character(len=*), parameter :: testfile = "test_output.dat"
        real(8), dimension(3) :: vec = [1.0d0, 2.0d0, 3.0d0]
        logical :: exists
        
        call write_data(testfile, vec)
        inquire(file=testfile, exist=exists)
        
        if (exists) then
            print *, "PASS: File was created."
            call execute_command_line("rm -f "//testfile)
        else
            print *, "FAIL: File was not created."
            stop 1
        end if
    end subroutine test_write_data_creates_file

end program test_mod_data_io
