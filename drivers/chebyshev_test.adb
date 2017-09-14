with Chebyshev, Ada.Text_IO, Numerics;
use  Chebyshev, Ada.Text_IO, Numerics;

procedure Chebyshev_Test is
   use Real_IO, Int_IO, Real_Functions;
   N : constant Nat  := 60;
   L : constant Real := -1.0e-2;
   R : constant Real := 1.1e-2;
   X : Real_Array := Chebyshev_Gauss_Lobatto (N, L, R);
   F : Real_Array (X'Range);
   A : Real_Array (X'Range);
   Y : Real_Array (1 .. 16);
   Z : Real;
   W : Real;
   
   function Test (X : in Real) return Real is
   begin
      --  return X * Exp (-X) * Cos (4.0 * X);
      return Sqrt (abs (X));
   end Test;
begin
   
   for I in X'Range loop
      F (I) := Test (X (I));
   end loop;

   A := CGL_Transform (F);
   
   for I in Y'Range loop
      Y (I) := L + (R - L) / Real (Y'Length - 1) * Real (I - 1);
   end loop;
   
   for I in Y'Range loop
      Z := Test (Y (I));
      W := Interpolate (A, Y (I), L, R);
      
      Put (Y (I), Aft => 3, Exp => 0); Put (",  ");
      Put (Z,     Aft => 5, Exp => 0); Put (",  ");
      Put (W,     Aft => 5, Exp => 0); Put (",  ");
      Put (100.0 * (Z - W) / (abs (Z) + 1.0e-10), Aft => 4, Exp => 0); New_Line;
   end loop;
   null;
end Chebyshev_Test;
