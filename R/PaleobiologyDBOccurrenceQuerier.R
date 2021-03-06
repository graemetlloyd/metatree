#' Palaeobiology Database Occurrence Querier
#'
#' @description
#'
#' Given a set of Paleobiology Database taxon name(s) or number(s) returns occurrence information for those tax(a).
#'
#' @param taxon_nos A vector of Paleobiology database taxon number(s) to retrieve from the database.
#' @param original Whether or not to return the original (TRUE) or resolved version (FALSE) of names.
#' @param breaker Size of breaker to use if querying a large number of taxa (reduces load on database of individual queries; default is 100).
#' @param RetainUncertainOccurrences Logical indicating whether or not to retain uncertain (i.e., aff, cf., ?, "") occurrences (defaults to FALSE).
#'
#' @details
#'
#' Uses the Paleobiology Database (\code{paleobiodb.org}) API (Peters and McLennen 2016) to query known taxon numbers and returns information on their occurrence as fossils (where this data is available in the database).
#'
#' @return
#'
#' A multi-column matrix with rows as occurrences.
#'
#' @author
#'
#' Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @references
#'
#' Peters, S. E. and McClennen, M., 2016. The Paleobiology Database application programming interface. \emph{Paleobiology}, \bold{42}, 1-7.
#'
#' @examples
#'
#' # Occurrence query for Allosaurus fragilis:
#' PaleobiologyDBOccurrenceQuerier(taxon_nos = "52962")
#'
#' @export PaleobiologyDBOccurrenceQuerier
PaleobiologyDBOccurrenceQuerier <- function(taxon_nos, original = TRUE, breaker = 100, RetainUncertainOccurrences = FALSE) {
  
  # SOME WAY TO ADD GEOGRAPHIC DATA FOR EXTANT TIPS?
  
  # Subfunction to break N numbers into breaker-sized blocks:
  NumberChunker <- function(N, breaker) {
    
    # Get total numebr of chunks:
    NChunks <- ceiling(N / breaker)
    
    # Make initial list of numbers:
    ListOfNumbers <- rep(list(1:breaker), NChunks)
    
    # Update last element of list (which can be less than breaker) into the correct size (if required):
    if((N %% breaker) > 0) ListOfNumbers[[length(ListOfNumbers)]] <- 1:(N %% breaker)
    
    # Add cumulative breaker value to all numbers (if required):
    if(N > breaker) ListOfNumbers <- mapply('+', ListOfNumbers, c(0, cumsum(rep(breaker, NChunks - 1))), SIMPLIFY = FALSE)
    
    # Return list of numbers:
    return(ListOfNumbers)
    
  }
  
  # Build list of numbers to query:
  NumbersToQuery <- lapply(NumberChunker(N = length(taxon_nos), breaker = breaker), function(x) taxon_nos[x])
  
  # Build HTTP string(s):
  ResolvedHTTPStrings <- lapply(NumbersToQuery, function(x) ifelse(original, paste("https://paleobiodb.org/data1.2/occs/list.json?taxon_id=", paste(paste("var:", x, sep = ""), collapse = ","), "&show=coords,paleoloc", sep = ""), paste("https://paleobiodb.org/data1.2/occs/list.json?taxon_id=", paste(paste("txn:", x, sep = ""), collapse = ","), "&show=coords,paleoloc", sep = "")))
  
  # Get resolved json strings for each chunk:
  ResolvedJSON <- lapply(ResolvedHTTPStrings, function(x) {
    
    # Set resolved json to NA (used later to check results are coming back from server):
    resolvedjson <- NA
    
    # Set start value for counter (used later to avoid infinite loop):
    counter <- 0
    
    # While server has not been reached (and querying a taxon number):
    while(is.na(resolvedjson[[1]][1])) {
      
      # Attempt to acquire resolved taxon string:
      try(resolvedjson <- readLines(x), silent = TRUE)
      
      # If server was not reached:
      if(is.na(resolvedjson[[1]][1])) {
        
        # Update counter to record how many attempts to reach server have been made:
        counter <- counter + 1
        
        # If repeatedly failing to get results stop trying:
        if(counter == 100) stop("Server not responding after 100 straight attempts")
        
        # Wait two seconds before next attempt (also avoids overloading server):
        Sys.sleep(2)
        
      }
      
    }
    
    # Return resolvedjson string:
    return(resolvedjson)
    
  })
  
  # Extract data from json data:
  Output <- do.call(rbind, lapply(ResolvedJSON, function(x) {
    
    # Subfunction to extract specific parameter from json string:
    ParameterExtraction <- function(jsonstring, parameterstring) {
      
      # Extract specific parameter:
      output <- unlist(lapply(as.list(jsonstring), function(x) ifelse(length(grep(parameterstring, x)) > 0, gsub("\"", "", strsplit(strsplit(x, parameterstring)[[1]][2], ",")[[1]][1]), NA)))
      
      # Return output:
      return(output)
      
    }
    
    # Find any taxon numbers not in database:
    UnknownTaxonHits <- grep("Unknown taxon", x)
      
    # If found stop and warn user:
    if(length(UnknownTaxonHits) > 0) stop(paste("The following taxon numbers were not found in the database: ", paste(unlist(lapply(strsplit(x[UnknownTaxonHits], split = "Unknown taxon '|'\""), '[', 2)), collapse = ", "), sep = ""))
    
    # Isolate record line:
    jsonstring <- x[(grep("\\[", x) + 1):(grep("\\]", x) - 1)]
    
    # List of parameters to extract (check API documentation for meaning):
    Parameters <- c("cid", "idn", "tna", "oei", "eag", "lag", "lng", "lat", "pln", "pla")
    
    # Compile output:
    output <- do.call(cbind, lapply(as.list(Parameters), function(x) gsub("col:", "", ParameterExtraction(jsonstring, parameterstring = paste("\"", x, "\":", sep = "")))))
    
    # Add names:
    colnames(output) <- c("CollectionNo", "IdentifiedName", "TaxonName", "Age", "MaxMa", "MinMa", "Longitude", "Latitude", "PalaeoLongitude", "PalaeoLatitude")
    
    # Return output:
    return(output)
    
  }))
  
  # If there are any unidentified
  if(any(is.na(Output[, "IdentifiedName"]))) {
    
    # Store rows with NAs for identified name(s):
    Rows <- which(is.na(Output[, "IdentifiedName"]))
    
    # Overwrite NAs with taxon name:
    Output[Rows, "IdentifiedName"] <- Output[Rows, "TaxonName"]
    
  }
  
  # If do not want to retain uncertain occurrences:
  if(!RetainUncertainOccurrences) {
    
    # Identify any uncertain occurrences:
    UncertainOccurrences <- grep("\\\\|cf\\.|aff\\.|\\?", Output[, "IdentifiedName"])
    
    # If found, remove these from the output:
    if(length(UncertainOccurrences) > 0) Output <- Output[-UncertainOccurrences, , drop = FALSE]
    
  }

  # Perform a taxon query (to check for extant and no occurrence taxa):
  TaxonQuery <- PaleobiologyDBDescendantFinder(taxon_nos = taxon_nos)
  
  # If extant taxa are found add them to the matrix:
  if(any(sort(TaxonQuery[, "Extant"] == "1"))) Output <- rbind(Output, do.call(rbind, lapply(as.list(TaxonQuery[TaxonQuery[, "Extant"] == "1", "TaxonName"]), function(x) c("0", x, x, "Extant", "0", "0", NA, NA, NA, NA))))
  
  # Find any taxa without (definite) occurrences (excludes extant taxa):
  NoOccurrenceTaxa <- setdiff(TaxonQuery[, "TaxonName"], unique(Output[, "TaxonName"]))
  
  # If no occurrence are found add them to the matrix with NAs:
  if(length(NoOccurrenceTaxa) > 0) Output <- rbind(Output, do.call(rbind, lapply(as.list(NoOccurrenceTaxa), function(x) c("0", x, x, NA, NA, NA, NA, NA, NA, NA))))

  # Return output:
  return(Output)
  
}

# Triassic dinosaur species:
#taxon_nos <- c("251932", "64409", "66344", "159719", "192929", "182713", "85754", "66909", "174879", "172073", "66661", "55001", "55000", "56379", "133339", "54097", "65490", "96666", "56382", "162588", "370732", "370734", "372124", "55477", "67921", "91403", "378703", "347522", "57454", "243275", "68124", "370986", "242722", "153776", "243373", "119230", "133334", "131589", "57406", "142534", "66571", "53365", "57413", "57412", "56627", "62959", "243277", "378960", "56592", "65100", "56683", "57418", "56596", "64684", "68584", "56585", "327194", "56614", "56586", "66682", "68588", "55644", "322708", "54980", "335603")

# Testudinata species:
#taxon_nos <- c("253273", "178076", "252755", "175826", "252659", "361368", "63142", "376586", "316508", "362616", "80840", "83314", "317850", "317845", "317853", "317849", "317851", "316873", "80841", "288218", "239568", "362607", "362606", "365789", "128353", "128352", "127991", "255124", "127985", "128346", "238132", "231906", "259783", "288227", "125121", "379435", "239570", "128038", "150512", "173417", "353086", "364099", "312334", "361541", "235088", "235086", "235083", "128037", "96952", "318682", "235269", "90355", "63881", "125124", "77320", "241951", "288216", "165653", "165652", "358066", "288214", "173399", "288961", "86857", "375082", "365636", "375083", "316520", "359362", "359359", "255080", "255081", "312335", "173402", "363169", "322199", "362853", "362244", "264912", "316510", "264913", "63978", "83307", "318120", "173404", "83305", "359400", "151651", "358387", "358384", "358386", "359401", "235428", "235427", "366281", "235426", "377282", "366282", "252453", "235259", "280081", "86275", "365809", "377274", "359391", "337573", "358385", "362856", "235092", "359389", "359386", "241953", "235267", "359397", "344463", "344462", "313940", "98332", "342729", "179131", "231911", "197866", "179138", "231912", "365634", "273630", "128700", "386802", "286869", "361531", "83033", "231900", "128348", "231908", "231919", "231940", "231901", "231907", "105425", "235813", "105427", "235665", "105426", "236083", "316506", "56476", "253836", "229121", "235672", "235673", "235669", "235668", "259494", "361411", "380803", "253376", "265932", "253269", "253267", "256142", "340791", "235646", "288123", "235643", "265929", "235645", "253268", "340791", "258618", "243386", "258620", "235808", "235807", "258632", "234783", "235810", "297501", "235812", "258628", "374431", "99472", "268066", "345816", "235448", "364012", "235454", "268052", "357838", "253053", "359403", "359353", "235663", "359352", "365896", "309527", "249577", "260255", "253407", "323691", "238155", "356323", "253403", "258623", "254130", "253593", "364012", "254123", "235458", "254042", "253598", "362193", "317945", "364011", "365811", "235485", "235536", "362144", "362138", "254030", "242706", "253386", "268055", "253602", "235828", "375804", "212674", "374099", "376553", "376551", "376552", "253263", "340873", "376554", "268028", "375807", "343442", "337694", "360823", "243533", "235123", "360704", "361359", "360695", "172987", "360681", "289017", "289460", "361350", "289015", "360669", "361344", "347392", "361341", "289016", "361342", "243529", "65120", "360679", "217623", "361352", "205690", "360706", "361357", "361343", "235119", "361348", "235132", "360708", "360670", "360716", "54272", "360712", "320931", "360715", "361354", "361349", "361355", "360705", "361356", "205688", "289014", "361347", "54273", "352002", "376555", "360730", "376556", "352000", "377688", "377687", "375806", "289010", "351995", "268024", "268025", "351998", "351992", "351993", "351994", "364626", "364631", "235114", "121524", "172994", "121522", "361870", "121526", "243545", "375808", "360766", "131562", "358273", "119641", "131565", "131563", "288223", "352420", "67642", "374265", "373846", "254324", "254325", "253384", "316768", "66832", "290854", "65858", "296183", "184480", "235541", "234884", "63644", "253701", "57643", "63643", "360569", "212719", "67760", "203067", "212675", "67312", "212681", "342862", "202999", "259775", "344452", "288843", "66154", "67897", "112785", "67759", "238233", "65789", "67590", "238227", "241944", "67762", "56158", "213506", "342861", "347364", "67761", "289458", "82832", "358608", "153036", "360566", "253374", "253369", "253368", "239899", "361414", "288819", "358274", "235832", "202717", "239574", "307632", "241956", "243313", "67765", "367942", "288973", "68208", "235823", "120459", "57306", "376426", "238551", "235761", "235763", "317883", "253591", "374126", "235834", "65382", "235836", "197700", "340818", "258839", "65906", "197699", "238723", "202723", "136765", "259792", "65385", "65386", "377131", "331283", "81007", "238472", "248475", "197851", "345928", "235652", "234896", "234794", "286867", "156175", "82838", "235776", "361416", "105819", "82836", "345921", "243318", "345922", "243333", "82835", "54274", "235794", "235774", "279458", "172490", "313594", "243331", "179128", "173225", "230534", "256143", "235792", "235795", "82837", "235772", "243320", "387333", "235781", "235780", "235779", "197948", "202903", "202901", "156173", "235785", "235662", "172719", "343386", "172715", "276090", "172716", "360992", "377464", "377460", "377461", "360997", "376496", "377470", "377463", "360993", "377467", "377466", "105336", "377462", "118952", "377465", "377468", "376497", "376493", "360996", "360994", "53091", "376494", "360995", "376495", "377469", "377445", "377443", "377444", "361001", "361000", "376492", "172718", "105391", "172721", "172720", "337692", "172849", "172862", "360984", "243342", "243326", "172866", "65794", "172857", "172861", "234957", "360522", "322258", "172489", "96555", "339105", "238112", "243346", "373841", "374596", "360534", "347000", "374588", "346778", "360540", "339106", "346777", "364829", "360514", "236085", "247616", "317855", "351610", "351990", "316913", "358158", "351988", "351986", "351985", "351987", "238462", "351983", "264152", "264155", "322911", "238460", "264162", "238456", "259499", "92235", "287375", "92236", "361340", "253839", "253840", "361339", "361338", "308281", "259497", "259496", "259495", "361390", "259501", "67667", "128338", "127983", "235255", "235257", "334106", "365790", "235149", "235096", "235094", "127999", "127996", "128355", "128357", "173397", "368830", "128003", "128358", "231892", "213913", "165647", "364499", "128350", "241955", "263699", "244998", "289153", "365779", "83311", "105049", "247704", "291353", "291354", "291355", "172868", "358115", "358117", "377437", "377436", "105357", "113507", "358113", "331593", "358156", "339851", "358108", "358112", "253401", "364317", "231916", "319457", "172850", "319459", "238603", "319461", "366280", "132780", "307884", "341553", "116706", "344584", "344587", "132761", "375519", "96381", "366287", "121688", "173138", "121679", "366294", "366288", "366299", "366303", "364319", "366300", "347876", "358808", "347900", "347951", "380349", "358120", "235421", "235423", "235425", "374260", "374267", "373909", "373911", "377705", "374264", "374271", "203535", "376560", "375454", "364854", "376559", "365357", "382656", "261149", "341249", "261150", "377727", "377726", "377728", "377725", "377729", "105280", "105292", "105278", "105294", "105296", "105290", "366305", "235304", "377721", "132800", "374270", "375442", "373907", "375444", "376119", "203530", "376120", "132777", "203526", "238691", "376122", "132797", "376117", "376118", "132798", "376121", "376123", "238692", "366283", "341253", "341245", "341252", "346199", "347872", "376557", "365783", "375445", "376558", "377720", "338366", "338751", "253397", "377482", "377483", "376561", "203532", "347700", "341554", "105213", "105238", "366284", "105217", "366285", "105233", "105221", "366286", "105227", "105212", "105219", "364826", "86855", "235282", "235280", "374266", "382658", "341265", "376845", "341263", "235272", "364325", "359388", "235278", "346198", "341270", "376852", "377677", "238489", "364835", "337848", "347893", "364508", "341266", "376584", "235275", "380346", "96498", "341261", "364836", "364472", "172497", "364329", "377678", "238488", "364328", "341259", "365833", "383995", "365639", "235273", "235274", "383996", "347955", "364837", "341260", "376850", "377676", "341267", "364327", "376853", "346352", "341271", "365832", "341258", "341268", "235277", "374257", "364839", "364326", "376847", "172496", "341264", "363763", "346197", "364502", "235276", "364838", "358990", "365780", "341244", "347697", "377484", "374258", "373908", "235299", "373912", "235145", "235091", "235144", "359402", "179145", "234797", "358124", "346130", "373906", "373905", "377486", "364331", "364323", "377493", "377491", "358438", "374950", "377489", "356116", "242703", "377494", "242702", "242704", "351851", "377490", "242701", "364321", "358470", "377488", "364322", "377487", "377492", "377485", "364957", "344586", "338368", "173161", "339828", "338918", "339556", "339560", "339855", "339845", "339570", "339571", "340451", "340443", "377707", "358360", "358365", "358363", "358364", "339104", "377706", "358367", "68432", "316516", "373847", "365748", "365781", "365749", "377714", "365777", "377710", "377711", "377712", "377708", "377713", "377709", "340399", "135030", "238487", "172491", "380345", "58731", "332929", "358129", "358130", "235295", "331594", "350831", "235296", "373879", "375452", "373876", "373877", "377441", "373874", "373880", "373875", "377440", "373878", "341256", "377442", "341257", "373873", "377683", "377680", "377684", "377682", "377681", "337289", "380344", "105390", "96552", "366271", "366269", "105299", "96553", "366270", "173236", "366272", "375634", "377691", "377692", "373870", "316514", "340709", "339843", "340474", "377052", "340696", "137548", "365776", "137557", "340712", "316300", "380313", "340698", "137553", "340471", "377693", "276089", "377697", "377701", "377698", "377696", "377695", "377700", "377699", "377694", "377690", "377689", "235407", "235401", "381655", "339576", "339764", "236086", "235913", "236087", "341550", "65903", "65904", "65902", "57304", "217279", "217278", "197957", "261025", "235830", "235825", "261151", "238126", "238153", "238122", "57303", "238123", "253655", "339835", "340458", "340456", "338850", "339826", "338849", "363917", "338753", "316303", "339813", "339782", "373915", "373914", "375819", "375818", "341048", "365746", "373903", "373901", "373904", "373902", "105194", "96475", "96481", "361413", "377391", "377398", "121576", "121578", "377392", "366029", "365627", "377402", "377403", "377399", "377400", "377394", "121572", "121574", "377401", "173007", "121580", "366268", "377397", "235418", "316301", "96386", "364828", "364827", "119610", "90483", "90478", "90477", "119609", "100141", "100124", "364841", "90479", "100176", "229117", "383080", "173009", "105350", "339849", "366263", "375814", "375820", "207534", "121478", "105365", "384496", "105344", "365786", "374323", "366265", "366267", "377703", "377704", "53092", "121570", "105346", "364662", "377702", "53093", "105348", "383081", "172492", "380314", "235409", "376564", "376563", "376562", "341043", "376585", "375084", "341254", "172493", "373872", "375822", "80527", "235142", "374533", "375475", "340714", "105184", "143986", "373850", "366264", "375816", "105368", "137549", "105179", "378733", "338853", "374582", "364842", "364564", "319463", "186113", "81428", "361281", "110737", "364027", "361298", "375075", "364026", "358411", "375081", "361305", "361286", "373851", "358388", "375074", "234799", "347379", "56736", "203002", "377429", "377425", "377422", "377430", "377428", "377426", "377427", "377732", "377737", "377736", "377731", "377423", "377421", "377740", "377734", "377424", "377439", "336683", "373855", "373857", "373852", "235472", "373853", "373860", "361288", "373856", "373861", "361287", "373859", "373854", "373858", "361294", "361292", "373845", "373844", "361293", "373842", "356487", "377363", "377360", "377365", "377367", "375072", "361299", "375069", "361284", "361283", "375071", "375070", "375079", "375077", "375078", "375080", "235476", "375076", "375073", "81424", "81426", "238473", "191805", "81430", "387882", "365523", "344490", "131531", "131456", "238475", "131454", "344525", "344531", "255002", "255007", "255006", "255001", "255005", "365525", "342993", "344520", "239573", "344516", "268033", "322630", "346165", "345371", "343101", "343100", "346163", "343077", "337346", "366273", "343096", "344529", "131544", "365002", "179052", "143984", "344518", "172422", "363957", "344495", "344491", "344498", "358830", "344503", "358828", "363802", "344532", "344504", "92213", "344502", "344501", "373849", "344505", "358829", "364855", "317283", "128104", "373848", "253845", "344489", "344507", "344521", "86279", "361160", "356716", "356717", "356715", "344500", "131522", "131470", "131518", "131510", "301181", "131514", "131493", "131494", "131478", "131476", "344989", "239248", "98659", "371054", "355640", "371053", "131485", "131487", "131488", "131491", "131489", "131482", "243865", "238148", "388136", "197713", "375096", "375094", "197722", "131548", "374386", "179140", "131508", "179139", "131506", "131500", "131512", "131463", "153867", "131524", "131479", "131483", "131516", "131461", "131472", "131471", "348492", "363950", "131442", "131539", "211899", "131451", "131449", "131445", "131525", "131542", "373240", "344499", "238479", "247700", "351867", "351865", "217404", "234893", "234894", "55587", "65856", "235533", "332081", "352019", "57030", "352018", "253379", "352111", "253379", "352017", "261782", "270421", "235568", "113307", "247702", "253359", "347365", "238120", "174446", "243381", "253264", "168584", "137153", "351832", "131561", "131561", "344458", "253378", "351831", "180695", "178245", "253355")

#x <- PaleobiologyDBOccurrenceQuerier(taxon_nos)
#x <- matrix(cbind(rep(mean((as.numeric(x[, "MaxMa"]) + as.numeric(x[, "MinMa"])) / 2), nrow(x)), as.numeric(x[, "PalaeoLongitude"]), as.numeric(x[, "PalaeoLatitude"])), ncol = 3, dimnames = list(c(), c("recon_age", "paleolng", "paleolat")))
#maps <- getmap(ma = sort(unique(x[, "recon_age"])), model = "PALEOMAP", do.plot = FALSE)
#mapast(model = "PALEOMAP", data = as.data.frame(x), map = maps)
