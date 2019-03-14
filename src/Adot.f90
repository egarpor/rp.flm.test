      subroutine adot(n, inprod, Adot_vec)
      implicit none

      ! Arguments
      double precision inprod((n * n + n) / 2), Adot_vec((n * n - n + 2) / 2)
      integer n

      ! Local variables
      double precision aux, den, quo, sumr
      integer i, j, r, ij, ir, ii, jr, jj, rj, rr, aux_i, aux_j, aux_r, ind_ij
      real, parameter :: pi = acos(-1.0)

      ! The first element of Adot_vec is the common diagonal element
      Adot_vec(1) = pi * (n+1)

      ! The rest of the elements are the lower triangle matrix of Adot
      do i = 2, n
        do j = 1, i-1

          ! Do the sum on the r index
          sumr = 0
          do r = 1, n

            ! From the definition of Aijr0
            if((i == r) .or. (j == r)) then

              ! Sum variable
              sumr = sumr + pi

            else

              ! Auxiliar variables for the indexes
              aux_i = i * (i - 1) / 2
              aux_j = j * (j - 1) / 2
              aux_r = r * (r - 1) / 2

              ! Indexes
              ij = aux_i + j
              ii = aux_i + i
              jj = aux_j + j
              rr = aux_r + r

              ! Lower triangular part of the matrix
              if (i > r) then

                ir = aux_i + r

              ! Upper triangular part
              else

                ir = aux_r + i

              end if

              ! Lower triangular part of the matrix
              if (r > j) then

                rj = aux_r + j

              ! Upper triangular part
              else

                rj = aux_j + r

              end if

              ! Symmetry
              jr = rj

              ! Computation of the quotient
              aux = inprod(ij) + inprod(rr)
              den = sqrt((aux - 2 * inprod(ir)) * (aux - 2 * inprod(jr)))
              quo = (aux - inprod(ir) - inprod(rj)) / den

              ! Avoid numerical problems on acos
              quo = max(quo, -1.0)
              quo = min(quo, 1.0)

            ! Sum
            sumr = sumr + abs(pi - acos(quo))

            end if

          end do

          ! Enter the ij-th element of A
          ind_ij = 1 + ((i - 1) * (i - 2) / 2) + j
          Adot_vec(ind_ij) = sumr

        end do
      end do

      return
      end
