with Numerics, Ada.Text_IO, Molecular_Dynamics;
use  Numerics, Ada.Text_IO, Molecular_Dynamics;

procedure Bubble is
   use Real_Functions, Real_IO, Int_IO;
   N  : constant Nat := 3;
   M  : constant Nat := 5;
   Dt : constant Real := 1.0e-4;
   U  : constant Pos2D := (0.2, 0.0);
   V  : constant Pos2D_Vector (1 .. M + 1) := (others => Dt * U);
   T_Final : constant Real := 5.0;
   T_Output : constant Real := 0.1;
   
   File : File_Type;
   
   R : Pos2D_Vector := Initialize_Lattice (N, M, 4 * (M + 1) + 2);
   S, R_New : Pos2D_Vector := R;
   T : Real := 2.0 * Dt;
   Tmp : Real;
   Iter : Nat := 1;
begin
   
   
   Create (File => File, Name => "out.xyz");
   
   -- Initializations ----------
   Output (R, File);
   R_New (1 .. M + 1) := R_New (1 .. M + 1) + V (1 .. M + 1);
   Output (R_New, File);
   -----------------------------
   while (T < T_Output) loop
      S := Verlet (R => R,
		   R_New => R_New,
		   Dt => Dt);
      R     := R_New;
      R_New := S;
      T     := T + Dt;
   end loop;
   Output (S, File);
   
   while (T < T_Final) loop
      Put ("Iter = "); Put (Iter); New_Line; Iter := Iter + 1;
      Put ("Time: "); Put (T); New_Line;
      Tmp := T + T_Output;
      while (T < Tmp) loop
	 S := Verlet (R => R,
		      R_New => R_New,
		      Dt => Dt);
	 R     := R_New;
	 R_New := S;
	 T     := T + Dt;
      end loop;
      Output (S, File);
   end loop;
   
   
   
   Close (File => File);
   
end Bubble;
