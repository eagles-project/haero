module validation

contains
  ! Given the name of a Skywalker input YAML file, determine the name of the
  ! corresponding output Python module file.
  ! @param [in] input_file The name of the Skywalker input YAML file.
  function output_name(input_file) result(filename)
    character(len=*), intent(in) :: input_file

    integer :: suffix, slash
    character(len=255) :: filename

    ! Generate an output file name based on the name of the input file.
    suffix = index(trim(input_file), ".yaml")
    if (suffix > 0) then
      ! Look for the last slash in the filename.
      slash = index(trim(input_file), "/", back=.true.)
      if (slash == 0) slash = 1
      filename = "haerof90_" // input_file(slash+1:suffix-1) // ".py"
    else
      filename = trim(input_file)
    end if
  end function

end module
