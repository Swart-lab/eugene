/**
   This is the menu creation code - place it right after you body tag
   Feel free to add this to a stand-alone js file and link it to your page.
**/

//Menu object creation
oCMenu=new makeCM("oCMenu") //Making the menu object. Argument: menuname

//Menu properties   
oCMenu.pxBetween=17
oCMenu.fromLeft=5 
oCMenu.fromTop=160
oCMenu.rows=0 
oCMenu.menuPlacement=0

oCMenu.onlineRoot="http://www.inra.fr/bia/T/EuGene/" 
oCMenu.resizeCheck=0 
oCMenu.wait=1000 
oCMenu.fillImg=""
oCMenu.zIndex=0
      
//Background bar properties
oCMenu.useBar=1
oCMenu.barWidth="menu"
oCMenu.barHeight="menu" 
oCMenu.barClass="clBar"
oCMenu.barX="menu"
oCMenu.barY="menu"
oCMenu.barBorderX=0
oCMenu.barBorderY=0
oCMenu.barBorderClass=""

//Level properties - ALL properties have to be spesified in level 0
oCMenu.level[0]=new cm_makeLevel() //Add this for each new level
oCMenu.level[0].width=103
oCMenu.level[0].height=25
oCMenu.level[0].regClass="clLevel0"
oCMenu.level[0].overClass="clLevel0over"
oCMenu.level[0].borderX=1
oCMenu.level[0].borderY=1
oCMenu.level[0].borderClass="clLevel0border"
oCMenu.level[0].offsetX=0 
oCMenu.level[0].offsetY=0
oCMenu.level[0].cols=0
oCMenu.level[0].rows=0
oCMenu.level[0].align="right" 

//EXAMPLE SUB LEVEL[1] PROPERTIES - You have to spesify the properties you
//want different from LEVEL[0] - If you want all items to look the same just
//remove this
oCMenu.level[1]=new cm_makeLevel()
oCMenu.level[1].width=158
oCMenu.level[1].height=22
oCMenu.level[1].regClass="clLevel1"
oCMenu.level[1].overClass="clLevel1over"
oCMenu.level[1].style=""
oCMenu.level[1].align="right" 
oCMenu.level[1].offsetX=0
oCMenu.level[1].offsetY=0
oCMenu.level[1].borderClass="clLevel1border"
oCMenu.level[1].borderX=1 
oCMenu.level[1].borderY=1
oCMenu.level[1].rows=0
oCMenu.level[1].align="right" 

/******************************************
      Menu item creation:
      myCoolMenu.makeMenu(name, parent_name, text, link, target, width, height,
                          regImage, overImage, regClass, overClass , align,
                          rows, nolink, onclick, onmouseover, onmouseout) 
*************************************/
oCMenu.makeMenu('top0','','<img src="Images/m_home.jpg">','index.html')

oCMenu.makeMenu('top1','','<img src="Images/m_plugins.jpg">','plugins.html')
  
oCMenu.makeMenu('top2','','<img src="Images/m_orga.jpg">','organisms.html')

oCMenu.makeMenu('top2b','','<img src="Images/m_references.jpg">','references.html')
    
oCMenu.makeMenu('top3','','<img src="Images/m_doc.jpg">','https://mulcyber.toulouse.inra.fr/docman/view.php/10/1047/EuGeneDoc.pdf')

oCMenu.makeMenu('top4','','<img src="Images/m_download.jpg">','https://mulcyber.toulouse.inra.fr/frs/?group_id=10')

oCMenu.makeMenu('top5','','<img src="Images/m_gforge.jpg">','https://mulcyber.toulouse.inra.fr/projects/eugene/')
   
oCMenu.makeMenu('top6','','<img src="Images/m_links_contacts.jpg">','links.html')
  
//Leave this line - it constructs the menu
oCMenu.construct()
