with Numerics, Ada.Text_IO;
use  Numerics, Ada.Text_IO;

procedure Landscape is
   use Real_IO, Int_IO, Real_Functions;
   α : Real := 1.0e2;
   
   function Phi (R : in Real) return Real is
   begin
      return 0.5 * (1.0 + Tanh (7.0 * (R - 0.8)));
   end Phi;
   
   function PE (Q : in Real_Vector) return Real is
      PE_G, PE_M : Real;
      T     : Real renames Q (1);
      S     : Real renames Q (2);
      R     : constant Real := 2.0 * abs (Cos (0.5 * (T + S)));
      Cs    : constant Real := Cos (S);
      Ct    : constant Real := Cos (T);
      C2tps : constant Real := Cos (2.0 * T +       S);
      Ctp2s : constant Real := Cos (T       + 2.0 * S);
      Vo, Vi, Tmp, Ro : Real;
   begin
      Ro := 0.7;
      Tmp := 1.0 / R;
      PE_G := Ct + 2.0 * C2tps;
      --  Vi := 1.0 / (R - Ro) ** 2;
      Vi := 0.7 * Exp (-0.5 * (R - Ro)) / (R - Ro);
      Vo := (Tmp ** 3) * Cos (2.0 * (T + S))
	- 3.0 * (Tmp ** 5) * ((Ct + C2tps) * (Cs + Ctp2s));
      PE_M := Cos (2.0 * T) - 3.0 * Ct ** 2
	+  Cos (2.0 * S) - 3.0 * Cs ** 2 + Vo * Phi (R) + Vi * (1.0 - Phi (R));
      return ((α / 6.0) * PE_M + PE_G);
   end PE;
   R : constant Real := π * 0.73; -- * 2.0 / 3.0;
   
   function PE2 (Q : in Real_Vector) return Real is
      PE_G, PE_M : Real;
      T     : Real renames Q (1);
      S     : Real renames Q (2);
      R     : constant Real := 2.0 * abs (Cos (0.5 * (T + S)));
      Cs    : constant Real := Cos (S);
      Ct    : constant Real := Cos (T);
      C2tps : constant Real := Cos (2.0 * T +       S);
      Ctp2s : constant Real := Cos (T       + 2.0 * S);
      Vo, Tmp : Real;
   begin
      Tmp := 1.0 / R;
      PE_G := Ct + 2.0 * C2tps;
      Vo := (Tmp ** 3) * Cos (2.0 * (T + S))
	- 3.0 * (Tmp ** 5) * ((Ct + C2tps) * (Cs + Ctp2s));
      PE_M := Cos (2.0 * T) - 3.0 * Ct ** 2 + Cos (2.0 * S) - 3.0 * Cs ** 2 + Vo;
      return ((α / 6.0) * PE_M + PE_G);
   end PE2;
   
   
   function Coordinate_Transform (X : in Real_Vector) return Real_Vector is
      Y   : Real_Vector (1 .. 2);
   begin
      Y (1) := 0.5 *  (3.0 * X (1) + X (2));
      Y (2) := 0.5 * (-3.0 * X (1) + X (2));
      return Y;
   end Coordinate_Transform;
   
   
   N : constant Nat := 200;
   Dx : constant Real := 2.0 / Real (N);
   X, Y : Real_Vector (1 .. 2);
   Num : Pos;
   S : constant Real := 180.0 / π;
   Tmp : Real;
   File : File_Type;
begin
   
   Create (File, Name => "landscape.vtk");
   -- headers
   Put_Line (File, "# vtk DataFile Version 2.0");
   Put_Line (File, "landscape");
   Put_Line (File, "ASCII");
   New_Line (File);
   -- grid points (the corners)
   Num := (N + 1) * (N + 1);
   Put_Line (File, "DATASET UNSTRUCTURED_GRID");
   Put (File, "POINTS "); Put (File, Num, 0); Put_Line (File, " float");
   for I in 1 .. N + 1 loop
      X (1) := -1.0 + Real (I - 1) * Dx;
      for J in 1 .. N + 1 loop
	 X (2) := -1.0 + Real (J - 1) * Dx;
	 Y := R * Coordinate_Transform (X);
	 Put (File, S * Y (1)); Put (File, " ");
	 Put (File, S * Y (2));
   	 Y := R * Coordinate_Transform (X);
   	 Tmp := PE (Y);
   	 Put (File, " "); Put (File, Tmp); New_Line (File);
	 --  Put_Line (File, " 0.0");
      end loop;
   end loop;
   New_Line (File);
   -- identify cells
   Num := N * N; -- number of cells
   Put (File, "CELLS "); Put (File, Num, Width => 0); Put (File, " ");
   Num := 5 * Num; -- number of points total used to identify cells
   Put (File, Num, Width => 0); Put (File, " ");
   New_Line (File);
   for I in 1 .. N loop
      for J in 1 .. N loop
	 Put (File, "4 ");
	 Num := (I - 1) * (N + 1) + J - 1; -- bottom-left corner
	 Put (File, Num, Width => 0); Put (File, " ");
	 Num := Num + 1; -- bottom-right corner
	 Put (File, Num, Width => 0); Put (File, " ");
	 Num := I * (N + 1) + J; -- top-right corner
	 Put (File, Num, Width => 0); Put (File, " ");
	 Num := Num - 1; -- top-left corner
	 Put (File, Num, Width => 0); New_Line (File);
      end loop;
   end loop;
   New_Line (File);
   --- identify cell types
   Num := N * N; -- number of cells
   Put (File, "CELL_TYPES "); Put (File, Num, Width => 0); New_Line (File);
   for I in 1 .. Num loop Put_Line (File, "9"); end loop; -- 9 <--> quad
   New_Line (File);
   
   Num := N * N;
   Put (File, "CELL_DATA "); Put (File, Num, Width => 0); New_Line (File);
   Put_Line (File, "SCALARS PE_difference float 1");
   Put_Line (File, "LOOKUP_TABLE default");
   for I in 1 .. N loop
      X (1) := -1.0 + (Real (I) - 0.5) * Dx;
      for J in 1 .. N loop
   	 X (2) := -1.0 + (Real (J) - 0.5) * Dx;
   	 Y := R * Coordinate_Transform (X);
   	 Tmp := abs (PE (Y) - PE2 (Y));
   	 if 2.0 * abs (Cos (0.5 * (Y (1) + Y (2)))) < 1.0 then
   	    Put (File, 0.0); New_Line (File);
   	 else
   	    Put (File, Tmp); New_Line (File);
   	 end if;
      end loop;
   end loop;
   New_Line (File);
   
   --  Num := N * N;
   --  Put (File, "CELL_DATA "); Put (File, Num, Width => 0); New_Line (File);
   --  Put_Line (File, "SCALARS distance float 1");
   --  Put_Line (File, "LOOKUP_TABLE default");
   --  for I in 1 .. N loop
   --     X (1) := -1.0 + (Real (I) - 0.5) * Dx;
   --     for J in 1 .. N loop
   --  	 X (2) := -1.0 + (Real (J) - 0.5) * Dx;
   --  	 Y := R * Coordinate_Transform (X);
   --  	 Tmp := 2.0 * abs (Cos (0.5 * (Y (1) + Y (2))));
   --  	 Put (File, Tmp); New_Line (File);
   --     end loop;
   --  end loop;
   
   Close (File);
   
   
   --  Put_Line ("r, min_pe, max_pe");
   --  for I in 0 .. N - 1 loop
   --     R := Real (I) / Real (N - 1) * 0.999 * π;
   --     Max_Min (R, N, Max, Min);
   --     Put (R); Put (", ");
   --     Put (Min); Put (", ");
   --     Put (Max); New_Line;
   --  end loop;
   
end Landscape;
