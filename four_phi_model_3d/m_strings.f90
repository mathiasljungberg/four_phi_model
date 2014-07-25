module m_strings
contains
  
  subroutine string_int_concatenate(string, int)
    character(80), intent(inout):: string
    integer, intent(in):: int
    
    character(80):: string2
    
    write(string2,*) int
    string = trim(adjustl(string)) //  "_" // trim(adjustl(string2)) 
    
  end subroutine string_int_concatenate
  
  subroutine string_string_concatenate(string, string2)
    character(80), intent(inout):: string
    character(*), intent(in):: string2
    
    string = trim(adjustl(string)) // trim(adjustl(string2)) 
    
  end subroutine string_string_concatenate
  
  
end module m_strings
