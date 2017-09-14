with Chebyshev, Ada.Text_IO, Numerics;
use  Chebyshev, Ada.Text_IO, Numerics;

procedure Chebyshev_Test is
   use Real_IO, Int_IO, Real_Functions;
   N : constant Nat := 8;
   X : Real_Array := Chebyshev_Gauss_Lobatto (N, -1.0, 1.0);
   F : Real_Array (X'Range);
   A : Real_Array (X'Range);
   Y : Real_Array (1 .. 21);
   Z : Real;
begin
   
   for I in X'Range loop
      F (I) := (X (I)**2 - 0.5) * Sin (4.0 * X (I));
   end loop;

   A := CGL_Transform (F, -1.0, 1.0);
   
   for I in Y'Range loop
      Y (I) := -1.0 + 0.1 * Real (I - 1);
   end loop;
   
   for I in Y'Range loop
      Z := (Y (I) ** 2 - 0.5) * Sin (4.0 * Y (I));
      --  Put (I); Put (",  ");
      Put (Y (I), Aft => 3, Exp => 0); Put (",  ");
      Put (Z, Aft => 5, Exp => 0); Put (",  ");
      Put (Interpolate (A, Y (I)), Aft => 5, Exp => 0); Put (",  ");
      Put (Z - Interpolate (A, Y (I)), Aft => 5, Exp => 3); New_Line;
   end loop;
   null;
end Chebyshev_Test;
