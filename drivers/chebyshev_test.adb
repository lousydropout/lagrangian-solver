with Chebyshev, Ada.Text_IO, Numerics;
use  Chebyshev, Ada.Text_IO, Numerics;

procedure Chebyshev_Test is
   use Real_IO, Int_IO, Real_Functions;
   
   N : constant Nat  := 41;
   L : constant Real := -0.45;
   R : constant Real := 4.7;
   X : Real_Vector := Chebyshev_Gauss_Lobatto (N, L, R);
   F : Real_Vector (X'Range);
   A : Real_Vector (X'Range);
   Y : Real_Vector (1 .. 27);
   Z : Real;
   W : Real;

   function Test (X : in Real) return Real is
   begin
      return X * Exp (-X) * Cos (4.0 * X) + Sin (X) ** 2 / (X + 1.0);
   end Test;
   
begin
   
   
   
   
   for I in X'Range loop
      F (I) := Test (X (I));
   end loop;

   A := CGL_Transform (F);
   
   for I in Y'Range loop
      Y (I) := L + (R - L) * Real (I - 1) / Real (Y'Length - 1);
      Z     := Test (Y (I));
      W     := Interpolate (A, Y (I), L, R);
      
      Put (Y (I), Aft => 3, Exp => 0); Put (",  ");
      Put (Z,     Aft => 5, Exp => 0); Put (",  ");
      Put (W,     Aft => 5, Exp => 0); Put (",  ");
      Put (100.0 * (Z - W) / (abs (Z) + 1.0e-10), Aft => 4, Exp => 0); New_Line;
   end loop;
   
   null;
end Chebyshev_Test;
