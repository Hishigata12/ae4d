//Closes all open windows that are not images

while (null != ij3d.WindowManager.getCurrentImage()) { 
    img = ij3d.WindowManager.getCurrentImage(); 
    img.changes = false; 
    img.close(); 
}