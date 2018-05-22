//Closes all open windows that are not images
list = getList("window.titles"); 
     for (i=0; i<list.length; i++){ 
     winame = list[i]; 
     	selectWindow(winame); 
     run("Close"); 
}
//Closes the 3D viewer
call("ij3d.ImageJ3DViewer.close");