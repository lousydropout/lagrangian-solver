with Numerics;
use  Numerics;
package body Chebyshev is
   
   
   function CGL_Transform (F : in Real_Array;
			   L : in Real := -1.0;
			   R : in Real :=  1.0) return Real_Array is
      N : constant Nat := F'Length;
      X : constant Real_Array  := Chebyshev_Gauss_Lobatto (N);
      G : Real_Array := (2.0 / Real (N - 1)) * F ;
      T : Real_Matrix (1 .. N, 1 .. N);
   begin
      G (1) := 0.5 * G (1);
      G (N) := 0.5 * G (N);
      
      for J in X'Range loop
	 T (1, J) := 1.0; 
	 T (2, J) := X (J);
	 for I in 3 .. N loop
	    T (I, J) := 2.0 * X (J) * T (I - 1, J) - T (I - 2, J);
	 end loop;
      end loop;
      
      return (T * G);
   end CGL_Transform;
   
   function Chebyshev_Gauss_Lobatto (N : in Nat;
				     L : in Real := -1.0;
				     R : in Real :=  1.0) return Real_Array is
      use Real_Functions;
      K : constant Real := (R - L) / 2.0;
      X : Real_Array (1 .. N);
      Y : Real;
   begin
      for I in X'Range loop
	 Y     := Real (I - 1) * π / Real (N - 1);
	 X (I) := L + K * (1.0 - Cos (Y));
      end loop;
      return X;
   end Chebyshev_Gauss_Lobatto;
   
   
   function Derivative_Matrix (N    : in Nat;
			       L, R : in Real) return Real_Matrix is
      use Real_Functions;
      M : constant Pos         := N - 1;
      K : constant Real        := 2.0 / (R - L);
      X : constant Real_Array  := Chebyshev_Gauss_Lobatto (N);
      P : Real_Array  (1 .. N) := (others => 1.0);
      D : Real_Matrix (1 .. N, 1 .. N);
   begin
      P (1) := 2.0; P (N) := 2.0;
      -- Top-left and bottom-right corners
      D (1, 1) := -(1.0 + 2.0 * Real (M ** 2)) / 6.0;
      D (N, N) :=  (1.0 + 2.0 * Real (M ** 2)) / 6.0;
      -- Diagonals
      for I in D'First (1) + 1 .. D'Last (1) - 1 loop
	 D (I, I) := -0.5 * X (I) / (1.0 - X (I) ** 2);
      end loop;
      -- Non-diagonals      
      for I in D'Range (1) loop
	 for J in D'Range (2) loop
	    if I /= J then
	       D (I, J) := (P (I) / P (J)) / (X (J) - X (I));
	       if (I + J) mod 2 = 0 then D (I, J) := -D (I, J); end if;
	    end if;
	    D (I, J) := K * D (I, J); -- note: outside non-diag if-clause
	 end loop;
      end loop;
      
      return D;
   end Derivative_Matrix;

     
   function Interpolate (A : in Real_Array;
			 X : in Real) return Real is
      N : constant Nat := A'Length;
      T : Real_Array (1 .. N);
      Y : Real := 0.0;
   begin
      T (1) := 1.0; T (2) := X;
      for I in 3 .. N loop
	 T (I) := 2.0 * X * T (I - 1) - T (I - 2);
      end loop;
      T (1) := 0.5; T (N) := 0.5 * T (N);
      
      for I in 1 .. N loop
	 Y := Y + A (I) * T (I);
      end loop;
      return Y;
   end Interpolate;
end Chebyshev;
