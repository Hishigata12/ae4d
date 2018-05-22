import ij.*;
import ij.process.*;
import ij.gui.*;
import java.awt.*;
import ij.plugin.*;
import ij.plugin.frame.*;

public class get_Transform implements PlugIn {

	public void run(String arg) {
		final Transform3D t = new Transform3D();
		System.out.println("Hello");
		t.rotX(deg2rad(degree));

		//t.rotY(deg2rad(degree));

		//t.rotZ(deg2rad(degree));
	}

}
