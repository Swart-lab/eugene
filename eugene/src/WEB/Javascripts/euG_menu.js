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

oCMenu.offlineRoot="http:///www.inra.fr/bia/T/EuGene/" 
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
oCMenu.makeMenu('top0','','<img src="WEB/Images/m_home.jpg">','index.html')
oCMenu.makeMenu('sub00','top0','EuGène web site','index.html')
oCMenu.makeMenu('sub01','top0','Bia','http://www.inra.fr/bia/welcome.fr.shtml')
oCMenu.makeMenu('sub02','top0','Inra','http://www.inra.fr/')
oCMenu.makeMenu('sub03','top0','Génoplante','http://www.genoplante.com/htm/prehome.html')

oCMenu.makeMenu('top1','','<img src="WEB/Images/m_plugins.jpg">','plug_all.html')
oCMenu.makeMenu('sub10','top1','All plugins','plug_all.html')
oCMenu.makeMenu('sub11','top1','Signals','plug_sig.html')
  /*
    oCMenu.makeMenu('sub110','sub11','Start','plug_start.html')
    oCMenu.makeMenu('sub111','sub11','Splice','plug_splice.html')
    oCMenu.makeMenu('sub112','sub11','Stop','plug_stop.html')
    oCMenu.makeMenu('sub113','sub11','Transcript','plug_trans.html')
    oCMenu.makeMenu('sub114','sub11','FrameShift','plug_fs.html')
  */    
oCMenu.makeMenu('sub12','top1','Contents','plug_contents.html')
oCMenu.makeMenu('sub13','top1','Others','plug_others.html')
  
oCMenu.makeMenu('top2','','<img src="WEB/Images/m_orga.jpg">','orga_all.html')
oCMenu.makeMenu('sub20','top2','All organisms','orga_all.html')
oCMenu.makeMenu('sub21','top2','Arabidopsis thaliana','orga_ara.html')
oCMenu.makeMenu('sub22','top2','Oryza sativa','orga_ory.html')
oCMenu.makeMenu('sub23','top2','Medicago truncatula','orga_medi.html')
    
oCMenu.makeMenu('top3','','<img src="WEB/Images/m_doc.jpg">','doc.html')
oCMenu.makeMenu('sub30','top3','EuGène documentation','doc.html#top')
oCMenu.makeMenu('sub31','top3','EuGène screen shots','doc.html#screen')
oCMenu.makeMenu('sub32','top3','Bibliography','doc.html#bib')
oCMenu.makeMenu('sub33','top3','Slides','doc.html#slides')
  
oCMenu.makeMenu('top4','','<img src="WEB/Images/m_download.jpg">','http://carlit.toulouse.inra.fr/cgi-bin/EG')
oCMenu.makeMenu('top5','','<img src="WEB/Images/m_links_contacts.jpg">','links.html')
oCMenu.makeMenu('sub50','top5','Links','links.html#links')
oCMenu.makeMenu('sub51','top5','Contacts','links.html#contacts')
    
//Leave this line - it constructs the menu
oCMenu.construct()
