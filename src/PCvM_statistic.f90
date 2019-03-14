      subroutine pcvm_statistic(n, Adot_vec, residuals, statistic)
      implicit none

      ! Arguments
      double precision statistic, Adot_vec((n * n - n + 2) / 2), residuals(n)
      integer n

      ! Local variables
      double precision sums
      integer i, j, ind_ij

      ! Sum for the symmetric part
      sums = 0
      do i = 2, n
        do j = 1, i-1
          ind_ij = 1 + ((i - 1) * (i - 2) / 2) + j
          sums = sums + residuals(i) * Adot_vec(ind_ij) * residuals(j)
        end do
      end do

      ! Statistic computed as the sum of the diagonal and the symmetric part
      statistic = Adot_vec(1) * dot_product(residuals, residuals) + 2 * sums

      return
      end
