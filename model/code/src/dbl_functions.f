c***********************************************************************
c23456789012345678901234567890123456789012345678901234567890123456789012
c***********************************************************************
c   Double precision routine to convert a single precision variable
c   to a double precision variable (with trailing zeroes).
c***********************************************************************
c
        double precision function mydble(a)
        save
        real a
        character*30 chara

        write(chara,10)a
 10     format(e30.17)

        read(chara,10)mydble

        return
        end
c
c***********************************************************************
c***********************************************************************
c23456789012345678901234567890123456789012345678901234567890123456789012
c***********************************************************************
c   Double precision routine to convert a single precision
c   integer variable to a double precision floating variable
c   (with trailing zeroes).
c***********************************************************************
c
        double precision function dfloat(n)
        save
        integer n
        real a
        character*30 chara
c
        a = float(n)
        write(chara,10)a
 10     format(e30.17)

        read(chara,10)dfloat

        return
        end
c
c***********************************************************************

