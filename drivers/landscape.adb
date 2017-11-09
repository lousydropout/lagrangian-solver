with Numerics, Ada.Text_IO;
use  Numerics, Ada.Text_IO;

procedure Landscape is
   use Real_IO, Int_IO, Real_Functions;
   α : Real := 1.0e3;
   
   function Phi (R : in Real) return Real is
   begin
      return 0.5 * (1.0 + Tanh (50.0 * (R - 0.5)));
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
      Ro := 0.9;
      if R < Ro then return -250.0; end if;
      Tmp := 1.0 / R;
      PE_G := Ct + 2.0 * C2tps;
      Vi := 0.01 * Exp (-12.0 * (R - Ro)) / (R - Ro) ** 2;
      Vo := (Tmp ** 3) * Cos (2.0 * (T + S))
	- 3.0 * (Tmp ** 5) * ((Ct + C2tps) * (Cs + Ctp2s));
      PE_M := Cos (2.0 * T) - 3.0 * Ct ** 2
	+  Cos (2.0 * S) - 3.0 * Cs ** 2 
	+ Vo * Phi (R) + Vi * (1.0 - Phi (R));
      return ((α / 6.0) * PE_M + PE_G);
   end PE;
   R : constant Real := π * 0.7028;
   
   function PE2 (Q : in Real_Vector) return Real is
      PE_G, PE_M : Real;
      T     : Real renames Q (1);
      S     : Real renames Q (2);
      R     : constant Real := 2.0 * abs (Cos (0.5 * (T + S)));
      Cs    : constant Real := Cos (S);
      Ct    : constant Real := Cos (T);
      C2tps : constant Real := Cos (2.0 * T +       S);
      Ctp2s : constant Real := Cos (T       + 2.0 * S);
      Vo : Real;
   begin
      PE_G := Ct + 2.0 * C2tps;
      Vo := (Cos (2.0 * (T + S)) * R ** 2 - 3.0 * ((Ct + C2tps) * (Cs + Ctp2s)))
	/ R ** 5;
      PE_M := Cos (2.0 * T) - 3.0 * Ct ** 2 + Cos (2.0 * S) - 3.0 * Cs ** 2 + Vo;
      return ((α / 6.0) * PE_M + PE_G);
   end PE2;
   
   function PE3 (Q : in Real_Vector) return Real is
      PE_G, PE_M : Real := 0.0;
      T     : Real renames Q (1);
      S     : Real renames Q (2);
      Cs    : constant Real := Cos (S);
      Ct    : constant Real := Cos (T);
      C2tps : constant Real := Cos (2.0 * T +       S);
      Ctp2s : constant Real := Cos (T       + 2.0 * S);
      X, N : array (1 .. 3) of Real_Vector (1 .. 2);
      DX : Real_Vector (1 .. 2);
      Vo, R, Tmp : Real;
   begin
      N (1) := (0.0, 1.0);
      N (2) := (-Sin (2.0 * T), Cos (2.0 * T));
      N (3) := (-Sin (2.0 * (T + S)), Cos (2.0 * (T + S)));

      X (1) := (0.0, 0.0);
      X (2) := (-Sin (T), Cos (T));
      X (3) := (-Sin (T) - Sin (2.0 * T + S), 
		 Cos (T) + Cos (2.0 * T + S));
      
      R := 2.0 * abs (Cos (0.5 * (T + S)));
      Vo := Cos (2.0 * (T + S)) / R ** 3 
	- 3.0 * ((Ct + C2tps) * (Cs + Ctp2s)) / R ** 5;
      
      pragma Assert (abs (Norm (X (1) - X (2)) - 1.0) < 1.0e-15);
      pragma Assert (abs (Norm (X (3) - X (2)) - 1.0) < 1.0e-15);
      pragma Assert (abs (Norm (X (3) - X (1)) - R) < 1.0e-15);
      
      
      for I in 1 .. 3 loop
	 for J in I + 1 .. 3 loop
	    DX  := X (I) - X (J);
	    R   := Sqrt (Dot (DX, DX));
	    
	    Tmp := Dot (N (I), N (J)) * R ** 2 
	      - 3.0 * Dot (N (I), DX) * Dot (N (J), DX);
	    if I = 1 and J = 2 then
	       pragma Assert (abs (Tmp - (Cos (2.0 * T) - 3.0 * Ct ** 2))
				< 1.0e-15);
	    end if;
	    if I = 2 and J = 3 then
	       pragma Assert (abs (Tmp - (Cos (2.0 * S) - 3.0 * Cs ** 2)) 
				< 1.0e-15);
	    end if;
	    if I = 1 and J = 3 then
	       pragma Assert (abs (Tmp / R ** 5 - Vo) < 1.0e-15);
	    end if;
	    PE_M := PE_M + Tmp / R ** 5;
	 end loop;
      end loop;
      PE_G := Ct + 2.0 * C2tps;
      
      return ((α / 6.0) * PE_M + PE_G);
   end PE3;
   
   
   function Coordinate_Transform (X : in Real_Vector) return Real_Vector is
      Y   : Real_Vector (1 .. 2);
   begin
      Y (1) := 0.5 *  (3.0 * X (1) + X (2));
      Y (2) := 0.5 * (-3.0 * X (1) + X (2));
      return Y;
   end Coordinate_Transform;
   
   
   
   N  : constant Nat := 500;
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
   Put_Line (File, "SCALARS PE float 1");
   Put_Line (File, "LOOKUP_TABLE default");
   for I in 1 .. N loop
      X (1) := -1.0 + (Real (I) - 0.5) * Dx;
      for J in 1 .. N loop
   	 X (2) := -1.0 + (Real (J) - 0.5) * Dx;
   	 Y := R * Coordinate_Transform (X);
   	 --  Tmp := abs ((PE (Y) - PE2 (Y))); -- / PE2 (Y));
   	 --  if 2.0 * abs (Cos (0.5 * (Y (1) + Y (2)))) < 1.0 then
   	 --     Put (File, -1.0e-6); New_Line (File);
   	 --  else
   	 --     Put (File, Tmp); New_Line (File);
   	 --  end if;
	 Put (File, PE (Y)); New_Line (File);
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
   
end Landscape;
