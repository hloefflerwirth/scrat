pipeline.htmlSummary <- function(env)
{
  dirname <- paste(env$files.name, "- Results/html")
  dir.create(dirname, showWarnings=FALSE, recursive=TRUE) 
  
  filename <- file.path(paste(env$files.name, "- Results"), "html", "Summary.html")
  util.info("Writing:", filename)
  outfile <- file(filename, "w")
  
  cat("
    header { width: 100%; float: left; margin-bottom: 2em; }

    nav {  
      font-family: arial, sans-serif;
      position: relative;
      font-size:14px;
      color:#333;
      margin: 0px auto;
    }
    nav ul {
      list-style-type: none;
      border-bottom:4px solid #0056af;
      margin: 0;
    }
    nav ul li {
      float: left;
      position: relative;
      top:4px;
    }
    nav ul li a {
      text-align: center;
      padding:8px 12px;
      display:block;
      text-decoration:none;
      color:#333;
      background: #FFF;
    }
    nav ul li a.act {
      color:#FFF;
      background: #0056af;
    }
    nav ul li:hover > a {
      background: #0056af;
      color:#FFF;
    }
    nav ul li ul {
      display: none
    }
    nav ul li:hover ul {
      display: block;
      position: absolute;
      top:36px;
      padding:0;
      margin:0;
      width:auto;
      border:0;
    }
    nav ul li:hover ul li a {
      display:block;
      background:#FFF;
      color:#333;
      text-align: center;
      border: 1px solid #999;
      border-radius: 0;
      -moz-border-radius: 0;
      -webkit-border-radius: 0;
      font-weight:normal;
      width:120px;
    }
    nav ul li:hover ul li a:hover, nav ul li a:hover {
      background:#0056af;
      color:#FFF;
    }

    body { margin: 0; padding-left: 40px; padding-right: 40px; color: #333; background-color: #fff;
           font: normal normal normal 14px/180% sans-serif; }

    wrapper { width: 90%; min-width: 400px; max-width: 1000px; }
    h1, h2 { line-height: 210%; padding-top: 20px; border-bottom: 1px solid #eee; }
    dl { line-height: 180%; padding-top: 20px; }
    dl dt { width: 50%; float: left; color: #111; }
    dl dt:after { content: ':'; }

  ", sep="", file=outfile)

  close(outfile)
         
  #### Overview ####
  
  filename <- file.path(paste(env$files.name , "- Results"), "Summary.html")
  `outfile <- file(filename, "w")

  
  # <ul>
  #  <li><a target=\"_self\" href=\"Summary.html\">Samples</a></li>
  # </ul></li>

  cat("<html><head>
    <link rel=\"stylesheet\" type=\"text/css\" href=\"html/style.css\">
    <title>Summary of SCRAT analysis results - Overview</title>
    </head>

    <body>
      <header id=\"nav-menue\">
        <nav><ul>
          <li><a class=\"act\" target=\"_self\" href=\"Summary.html\" target=\"content\">Overview</a></li>
          <li><a target=\"_self\" href=\"html/safsari_block1.html\">I Preprocessing & quality check </a></li>
          <li><a target=\"_self\" href=\"html/safsari_block2.html\">II Group annotation & stability</a></li>
          <li><a target=\"_self\" href=\"html/safsari_block3.html\">III Portraying</a>
          <li><a target=\"_self\" href=\"html/safsari_block4.html\">IV Sample similarity analyses</a></li>
          <li><a target=\"_self\" href=\"html/safsari_block5.html\">V Pseudotime estimation</a></li>
        </ul></nav>
      </header>
      
      <div id=\"wrapper\">
        <h1>Overview</h1>
        <dl>
        <dt>Dataset name</dt>
        <dd>", env$preferences$dataset.name, "</dd>
        <dt>Number of samples</dt>
        <dd>", ncol(env$indata), "</dd>
        <dt>Number of groups</dt>
        <dd>", length(unique(env$group.labels)), "</dd>
        <dt>Number of genes</dt>
        <dd>", nrow(env$indata), "</dd>
        <dt>ID type of genes</dt>
        <dd>", ifelse( env$preferences$database.id.type!="", env$preferences$database.id.type, "not defined" ), "</dd>
        <dt>Dimension of the SOM</dt>
        <dd>", env$preferences$dim.1stLvlSom, " x ", env$preferences$dim.1stLvlSom, "</dd>
        <dt>Date</dt>
        <dd>", format(Sys.time(), "%a %b %d %X %Y %Z"), "</dd>
        <dt>Analyst</dt>
        <dd>", env$preferences$system.info["user"], "</dd>
        <dt>scrat version</dt>
        <dd>", env$preferences$session.info$otherPkgs$scrat$Version, "</dd>
        </dl>
      </div></body></html>", sep="", file=outfile)
  
  close(outfile)
  
 
  #### bock 1 ####
  
  filename <- file.path(paste(env$files.name , "- Results"), "html", "safsari_block1.html")
  outfile <- file(filename, "w")  
  
  cat("<!DOCTYPE html>
      <html>
      <head>
        <link rel=\"stylesheet\" type=\"text/css\" href=\"style.css\">
        <title>Summary of SCRAT analysis results - Preprocessing & quality check</title>
    	</head>

      <body>
    		<header id=\"nav-menue\">
          <nav><ul>
  					<li><a href=\"../Summary.html\">Overview</a></li>
  					<li><a class=\"act\" href=\"safsari_block1.html\">I Preprocessing & quality check</a></li>
  					<li><a href=\"safsari_block2.html\">II Group annotation & stability</a></li>
  					<li><a href=\"safsari_block3.html\">III Portraying</a></li>
  					<li><a href=\"safsari_block4.html\">IV Sample similarity analyses</a></li>
  					<li><a href=\"safsari_block5.html\">V Pseudotime estimation</a></li>
  				</ul></nav>
    		</header>

      <div id=\"wrapper\">
      <h1>Block I: Preprocessing & quality check</h1>
      <h2>Read count preprocessing</h2>", sep="", file=outfile)

  if(env$preferences$preprocessing$count.processing)
  {
    cat("Read counts are preprocessed using standard workflow:<br>
      Cells with library size 3 MAD below median are discarded, as well as cells with number of expressed genes 3 MAD blow median (see p.1 in PDF below).
      Read counts are then translated to expression values using scale factor normalization (p.2 shows the factors).

      <object data=\"../Data Overview/Read count preprocessing.pdf\" type=\"application/pdf\" style=\"width:100%;height:600px;border:1px solid black\"></object>", sep="", file=outfile)
  } else
  {
    cat("skipped", sep="", file=outfile)
  }

  cat("<h2>Input data overview</h2>
      Data distribution (p.1 in PDF below), mean expression vs. variance plot (p.2; inter-quartile ranges (IQR) are shown to identify potential outliers), and individual cells' and group-wise aggregated expression level (p.3) are provided.
      <object data=\"../Data Overview/Data Distribution.pdf\" type=\"application/pdf\" style=\"width:100%;height:600px;border:1px solid black\"></object>
 
      <h2>Cell cycle phase classification</h2>", sep="", file=outfile)

  if(file.exists(file.path(paste(env$files.name , "- Results"), "Data Overview", "Cell cycle phase.pdf")))
  {    
    cat("Cell cycle prediction was applied to classify cells into G1 (G1 score above 0.5), G2/M (G2/M score above 0.5) or S phase (otherwise) based on 
        expression patterns of pre-trained classifiers. 
        <object data=\"../Data Overview/Cell cycle phase.pdf\" type=\"application/pdf\" style=\"width:100%;height:600px;border:1px solid black\"></object>", sep="", file=outfile)
  } else
  {
    cat("skipped", sep="", file=outfile)
  }
  
  cat("</div></body></html>", sep="", file=outfile)
  
  close(outfile) 
  
  
  #### bock 2 ####
  
  filename <- file.path(paste(env$files.name , "- Results"), "html", "safsari_block2.html")
  outfile <- file(filename, "w")  
  
  cat("<!DOCTYPE html>
      <html>
      <head>
        <link rel=\"stylesheet\" type=\"text/css\" href=\"style.css\">
    		<title>Summary of SCRAT analysis results - Group annotation & stability</title>
    	</head>
    
    	<body>
    		<header id=\"nav-menue\">
          <nav><ul>
            <li><a href=\"../Summary.html\">Overview</a></li>
            <li><a href=\"safsari_block1.html\">I Preprocessing & quality check</a></li>
            <li><a class=\"act\" href=\"safsari_block2.html\">II Group annotation & stability</a></li>
            <li><a href=\"safsari_block3.html\">III Portraying</a></li>
            <li><a href=\"safsari_block4.html\">IV Sample similarity analyses</a></li>
            <li><a href=\"safsari_block5.html\">V Pseudotime estimation</a></li>
          </ul></nav>
        </header>
    
      <div id=\"wrapper\">
      <h1>Block II: Group annotation & stability</h1>
      <h2>Data driven group annotation</h2>", sep="", file=outfile)

  if(file.exists(file.path(paste(env$files.name , "- Results"), "Summary Sheets - Groups", "PAT-groups assignment.pdf")))
  {    
    cat("Groups are defined utilizing combinatorics of expression module activation (PATs). Optimum number of groups is determined using sum squared error (SSE) estimation for increasing group numbers (p.1 in PDF below): 
        when optimum is reached, increasing the group number does not entail a significant decrease in SSE. This \'elbow\' in SSE courve is automatically detected using parallel shift of the red line.<br>
        Numbers of cells assigned to the resulting groups are given as barplot (p.2).
        <object data=\"../Summary Sheets - Groups/PAT-groups assignment.pdf\" type=\"application/pdf\" style=\"width:100%;height:600px;border:1px solid black\"></object>", sep="", file=outfile)
  } else
  {
    cat("skipped", sep="", file=outfile)
  }

  cat("<h2>Group silhouette coefficients</h2>
      Sillhouette coefficient (S) represents correlation of a cell's expression state to the centroid of its group, in relation to correlation to the nearest of the centroids of the other groups (p.1 in PDF below; individual cell and group-wise aggregated S).
      Positive S represents compact groups, whereas negative S indicate overlapping groups.
      <object data=\"../Summary Sheets - Groups/Group assignment.pdf\" type=\"application/pdf\" style=\"width:100%;height:600px;border:1px solid black\"></object>", sep="", file=outfile)
  
  cat("</div></body></html>", sep="", file=outfile)
  
  close(outfile) 
  
  
  
 
  
  #### bock 3 ####
  
  filename <- file.path(paste(env$files.name , "- Results"), "html", "safsari_block3.html")
  outfile <- file(filename, "w")  
  
  cat("<!DOCTYPE html>
      <html>
      <head>
        <link rel=\"stylesheet\" type=\"text/css\" href=\"style.css\">
	    	<title>Summary of SCRAT analysis results - Portraying</title>
      </head>
      
      <body>
    		<header id=\"nav-menue\">
          <nav><ul>
            <li><a href=\"../Summary.html\">Overview</a></li>
            <li><a href=\"safsari_block1.html\">I Preprocessing & quality check</a></li>
            <li><a href=\"safsari_block2.html\">II Group annotation & stability</a></li>
            <li><a class=\"act\" href=\"safsari_block3.html\">III Portraying</a></li>
            <li><a href=\"safsari_block4.html\">IV Sample similarity analyses</a></li>
            <li><a href=\"safsari_block5.html\">V Pseudotime estimation</a></li>
          </ul></nav>
        </header>
      
      <div id=\"wrapper\">
      <h1>Block III: Portraying</h1>
      <h2>Single cell expression portraits</h2>", sep="", file=outfile)

      if(file.exists(file.path(paste(env$files.name , "- Results"),"Expression portraits.pdf")))
      {
        cat("Each cell's expression landscape is described by the metagene expression values. They are arranged according to the underlying SOM grid and visualized by an appropriate color gradient (red and blue colors indicate over- and underexpressed metagenes, respectively; green and intermediate colors indicate no or low differential expression).
            The color patterns emerge as smooth textures representing the fingerprint of transcriptional activity of each cell and are mutually comparable. 
            <object data=\"../Expression portraits.pdf\" type=\"application/pdf\" style=\"width:100%;height:600px;border:1px solid black\"></object>", sep="", file=outfile)
      } else
      {
        cat("skipped", sep="", file=outfile)
      }

      
  cat("<h2>Gene localization and information table</h2>
      This spreadsheet contains various information about the genes analysed, such as mean expression level, metagene position, chromosome and a short description.<br>
      <a href=\"../CSV Sheets/Gene localization.csv\" target=\"_blank\">
      Download spreadsheet (csv)</a>    
      
      <h2>Group portraits and reports</h2>
      Group-specific portraits are calculated and visualized as the mean expression value of each metagene averaged over all cells of the group.
      The report sheets additionally provide general information about differential expression and corresponding gene list, and enriched gene sets.
      <object data=\"../Summary Sheets - Groups/Report.pdf\" type=\"application/pdf\" style=\"width:100%;height:600px;border:1px solid black\"></object>

      <h2>Spot module reports</h2>
      Complementary approaches were applied to define expression modules using the metagene landscape. The reports comprise several visualizations and information of the modules obtained:<br>
      <ul> 
      <li>Overview about the landscape used to define modules (p.1, not available for k-means clusters)</li> 
      <li>Assignment of module labels and topmost enriched gene sets (p.2)</li>
      <li>Supervised heatmap of module expression (p.3)</li>
      <li>Module report sheets providing general information, expression profile across all cells, list of genes contained in the module, and enriched gene sets (p.4ff)</li> 
      </ul>
      Overexpression spot modules are defined by collecting all overexpressed metagenes in the individual cells and applying an expression threshold.
      <object data=\"../Summary Sheets - Modules/Overexpression Spots/Report.pdf\" type=\"application/pdf\" style=\"width:100%;height:600px;border:1px solid black\"></object>", sep="", file=outfile)
  
  if(file.exists(file.path(paste(env$files.name , "- Results"),"Summary Sheets - Modules","Group Overexpression Spots","Report.pdf")))
  {
    cat("<br><br>
        Group overexpression spot modules collect metagenes consistently overexpressed in the different groups, respectively.
        <object data=\"../Summary Sheets - Modules/Group Overexpression Spots/Report.pdf\" type=\"application/pdf\" style=\"width:100%;height:600px;border:1px solid black\"></object>", sep="", file=outfile)
  }
  
  cat("<br><br>
      k-means modules are defined by applying k-means clustering to the metagene expression data.
      <object data=\"../Summary Sheets - Modules/K-Means Cluster/Report.pdf\" type=\"application/pdf\" style=\"width:100%;height:600px;border:1px solid black\"></object>

      <br><br>
      Clusters in the distance map are defined by a threshold citerion. This map contains mean Euclidean distance of each metagene to its adjacent neighbors, and indicates clusters in data space projected to the SOM.
      <object data=\"../Summary Sheets - Modules/D-Cluster/Report.pdf\" type=\"application/pdf\" style=\"width:100%;height:600px;border:1px solid black\"></object>
      
    </div></body></html>", sep="", file=outfile)
  
  close(outfile) 
  
  
  
  
  #### bock 4 ####
  
  filename <- file.path(paste(env$files.name, "- Results"), "html", "safsari_block4.html")
  outfile <- file(filename, "w")  
  
  cat("<!DOCTYPE html>
      <html>
      
  	<head>
      <link rel=\"stylesheet\" type=\"text/css\" href=\"style.css\">
  		<title>Summary of SCRAT analysis results - Sample similarity analyses</title>
  	</head>
  
  	  <body>
    		<header id=\"nav-menue\">
          <nav><ul>
            <li><a href=\"../Summary.html\">Overview</a></li>
            <li><a href=\"safsari_block1.html\">I Preprocessing & quality check</a></li>
            <li><a href=\"safsari_block2.html\">II Group annotation & stability</a></li>
            <li><a href=\"safsari_block3.html\">III Portraying</a></li>
            <li><a class=\"act\" href=\"safsari_block4.html\">IV Sample similarity analyses</a></li>
            <li><a href=\"safsari_block5.html\">V Pseudotime estimation</a></li>
          </ul></nav>
        </header>
  
      <div id=\"wrapper\">
      <h1>Block IV: Sample similarity analyses</h1>
      Different sample similariy metrics and methods are applied to the dimension reduced expression module data and metadata, respectively. 
      Application of downstream analyses on aggregated data instead of single gene data was shown to increase representativeness and reduce noisiness.

      <h2>Supervised and clustered heatmaps</h2>
      Heatmaps of module (p.1-2) and metagene (p.3-4) expression data with supervised and hierarchically clustered sample ordering.
      <object data=\"../Sample Similarity Analysis/Hierarchical Clustering.pdf\" type=\"application/pdf\" style=\"width:100%;height:600px;border:1px solid black\"></object>
  
      <h2>Independent Component Analysis</h2>
      Independent component analysis (ICA) distributes cells along axes of most variability similar to principal component analysis. However, restriction to othogonal axes is omitted in ICA.
      ICA is applied to module (p.1-2) and metagene (p.3-4) expression data, where the first three components are shown in 3d and pairwise 2d scatterplots.
      <object data=\"../Sample Similarity Analysis/Component Analysis.pdf\" type=\"application/pdf\" style=\"width:100%;height:600px;border:1px solid black\"></object>
      
      <h2>Correlation Spanning Tree</h2>
      Correlation spanning tree represents module (p.1) and metagene (p.2) data as graph with cells as nodes connected to a spanning tree of maximal mutual correlation between connected nodes.
      <object data=\"../Sample Similarity Analysis/Correlation Spanning Tree.pdf\" type=\"application/pdf\" style=\"width:100%;height:600px;border:1px solid black\"></object>
  
      <h2>t-SNE</h2>
      t-distributed stochastic neighbor embedding (t-SNE) is a nonlinear dimensionality reduction technique projecting cells into a two-dimensional coordinate system. It is applied to module (p.1) and metagene (p.2) data.
      <object data=\"../Sample Similarity Analysis/t-SNE.pdf\" type=\"application/pdf\" style=\"width:100%;height:600px;border:1px solid black\"></object>
  
      </div></body></html>", sep="", file=outfile)
  
  close(outfile) 
  
  
  
  
  #### bock 5 ####
  
  filename <- file.path(paste(env$files.name, "- Results"), "html", "safsari_block5.html")
  outfile <- file(filename, "w")  
  
  cat("<!DOCTYPE html>
      <html>
   
    	<head>
        <link rel=\"stylesheet\" type=\"text/css\" href=\"style.css\">
    		<title>Summary of SCRAT analysis results - Pseudotime estimation</title>
    	</head>
    
    	<body>
    		<header id=\"nav-menue\">
          <nav><ul>
            <li><a href=\"../Summary.html\">Overview</a></li>
            <li><a href=\"safsari_block1.html\">I Preprocessing & quality check</a></li>
            <li><a href=\"safsari_block2.html\">II Group annotation & stability</a></li>
            <li><a href=\"safsari_block3.html\">III Portraying</a></li>
            <li><a href=\"safsari_block4.html\">IV Sample similarity analyses</a></li>
            <li><a class=\"act\" href=\"safsari_block5.html\">V Pseudotime estimation</a></li>
          </ul></nav>
        </header>


      <div id=\"wrapper\">
      <h1>Block V: Pseudotime estimation</h1>", sep="", file=outfile)
  
  if(file.exists(file.path(paste(env$files.name, "- Results"), "Pseudotime Analysis", "Trajectory report.pdf")))
  {    
    cat("Pseudotime is estimated using wanderlust algorithm. In brief, this approach generates an ensemble of K-nearest neighbor graphs based on the metagene expression landscape. For each graph, trajectories are calculated and subsequently averaged to obtain the final pseudotime trajectory.
      
      <h2>Pseudotime trajectory report</h2>
      The report presents information about pseudotime estimation results:<br>
      <ul> 
        <li>Correlation spanning tree and k-nearest neighbor graph utilized in wanderlust algorithm (p.1). Color code indicates pseudotime scores of the cells (see color key)</li> 
        <li>Metagene expression patterns as function of pseudotime (p.2): Significancy of pseudotime dependency is mapped to the SOM grid (topleft panel), and top-100 signigicant metagenes are masked (bottomleft) and their expression porfiles across the cells shown as heatmap visualization (right)</li>
        <li>Group specific pseudotime information (p.3): F-score and p-value for differential scores in the groups (topleft), pseudotime score ordered according to groups and to pseudotime trajectory (bottomleft), and mapping of the groups to k-nearest neighbor graph (right)</li>
        <li>Top-20 significantly pseudotime dependent genes (p.4ff): General information and localization in the SOM (topleft), differential expression ordered according to groups and to pseudotime trajectory (bottomleft), and mapping of expression values to k-nearest neighbor graph (right)</li>
        </ul>
      <object data=\"../Pseudotime Analysis/Trajectory report.pdf\" type=\"application/pdf\" style=\"width:100%;height:600px;border:1px solid black\"></object>
      
      <h2>Pseudotime score table</h2>
      This spreadsheet contains the list of cells ranked for increasing pseudotime scores along with group an PAT information.<br>
      <a href=\"../CSV Sheets/Pseudotime scores.csv\" target=\"_blank\">
      Download spreadsheet (csv)</a>    

      <h2>Pseudotime significance table</h2>
      This spreadsheet contains the list of genes ranked for decreasing pseudotime dependency (p-value and fdr in likelihood test of generalized additive model).<br>
        <a href=\"../CSV Sheets/Pseudotime trajectory genes.csv\" target=\"_blank\">
        Download spreadsheet (csv)</a> 

        ", sep="", file=outfile)
  } else
  {
    cat("skipped", sep="", file=outfile)
  }

  cat("
    </div>
  </body>
</html>", sep="", file=outfile)

  close(outfile)
}
