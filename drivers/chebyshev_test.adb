with Chebyshev, Ada.Text_IO, Numerics;
use  Chebyshev, Ada.Text_IO, Numerics;

procedure Chebyshev_Test is
   use Real_IO, Int_IO;
   X : Real_Array := Chebyshev_Gauss_Lobatto (3, -1.0, 1.0);
   D : Real_Matrix := Derivative_Matrix (3, -1.0, 1.0);
   F : Real_Array (X'Range);
   Y : Real_Array (X'Range);
   Z : Real_Array (X'Range) := (others => 0.0);
begin
   
   for I in X'Range loop
      F (I) := X (I);
   end loop;
   Y := D * F;
   
   for I in D'Range (1) loop
      for J in D'Range (2) loop
	 Put (D (I, J), Aft => 3, Exp => 0); Put (", ");
      end loop;
      New_Line;
   end loop;
   New_Line;
   
   for I in X'Range loop
      for J in X'Range loop
	 Z (I) := Z (I) + D (I, J) * F (J);
      end loop;
      Put (I); Put (",  ");
      Put (X (I), Aft => 3, Exp => 0); Put (",  ");
      Put (F (I), Aft => 3, Exp => 0); Put (",  ");
      Put (Y (I), Aft => 3, Exp => 0); Put (",  ");
      Put (Z (I), Aft => 3, Exp => 0); New_Line;
   end loop;
   null;
end Chebyshev_Test;
