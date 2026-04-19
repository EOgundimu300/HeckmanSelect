#' American Express Credit Card Data
#'
#' Data on 13,444 credit card applications received by American Express
#' in a single month in 1988.  Of the full sample, 10,499 applications
#' were approved and the subsequent 12 months of spending and default
#' behaviour were observed on approved applicants.  The dataset is used
#' in Greene (1998) to demonstrate sample-selection corrections in
#' credit scoring, and is a standard benchmark dataset for reject
#' inference in the credit-scoring literature.
#'
#' @format A data frame with 13,444 observations on 56 variables:
#' \describe{
#'   \item{CARDHLDR}{binary selection indicator: 1 if application was
#'         approved, 0 otherwise.}
#'   \item{DEFAULT}{binary outcome: 1 if cardholder defaulted, 0 otherwise.
#'         Observed only for approved applicants (\code{CARDHLDR = 1}).}
#'   \item{AGE}{age of applicant in years.}
#'   \item{ACADMOS}{months living at current address.}
#'   \item{ADEPCNT}{number of dependents.}
#'   \item{AEMPMOS}{months at current job.}
#'   \item{MAJORDRG}{number of major derogatory credit reports.}
#'   \item{MINORDRG}{number of minor derogatory credit reports.}
#'   \item{OWNRENT}{home ownership indicator (1 = owns, 0 = rents).}
#'   \item{APADMOS}{months since previous address.}
#'   \item{AMAMIND}{indicator for applicant's main income source.}
#'   \item{INCOME}{applicant's monthly income.}
#'   \item{SELFEMPL}{self-employment indicator.}
#'   \item{TRADACCT}{number of trade accounts.}
#'   \item{INCPER}{income per dependent.}
#'   \item{EXP_INC}{expenditure-to-income ratio.}
#'   \item{EXPRATIO}{alternative expenditure ratio.}
#'   \item{CPTOPNB}{number of bank-card trade lines.}
#'   \item{CPTOPNG}{number of general-purpose trade lines.}
#'   \item{CPT30C}{count of accounts 30-plus days delinquent.}
#'   \item{CPTF30}{fraction of accounts 30-plus days delinquent.}
#'   \item{CPTAVRV}{average revolving balance.}
#'   \item{CBURDEN}{credit burden indicator.}
#'   \item{BANKSAV}{savings-account indicator.}
#'   \item{BANKCH}{chequing-account indicator.}
#'   \item{BANKBOTH}{both savings and chequing indicator.}
#'   \item{CREDMAJR}{major credit-card indicator.}
#'   \item{CREDDEPT}{department-store credit indicator.}
#'   \item{CREDGAS}{gas-card credit indicator.}
#'   \item{ACBINQ}{number of recent credit inquiries.}
#'   \item{ACURTRD}{number of currently open trade lines.}
#'   \item{AADLINC}{additional-income indicator.}
#'   \item{BUYPOWER}{index of purchasing power.}
#'   \item{BUILDMTL}{building-materials spending.}
#'   \item{APPAREL}{apparel spending.}
#'   \item{AUTO}{automotive spending.}
#'   \item{DEPTSTOR}{department-store spending.}
#'   \item{DRUGSTOR}{drugstore spending.}
#'   \item{EATDRINK}{food-and-drink spending.}
#'   \item{FURN}{furniture spending.}
#'   \item{GAS}{gasoline spending.}
#'   \item{GROWTH}{income-growth indicator.}
#'   \item{CLERICAL}{clerical-occupation indicator.}
#'   \item{MGT}{management-occupation indicator.}
#'   \item{MILITARY}{military-service indicator.}
#'   \item{OTHERJOB}{other-occupation indicator.}
#'   \item{PROF}{professional-occupation indicator.}
#'   \item{SALES}{sales-occupation indicator.}
#'   \item{MEDAGE}{median age in applicant's ZIP code.}
#'   \item{MEDINC}{median income in applicant's ZIP code.}
#'   \item{PCTBLACK}{percent Black population in ZIP code.}
#'   \item{PCTSPAN}{percent Hispanic population in ZIP code.}
#'   \item{PCTCOLL}{percent college-educated in ZIP code.}
#'   \item{PCTEMPL}{percent employed in ZIP code.}
#'   \item{PCTOWN}{percent home-owners in ZIP code.}
#'   \item{UNEMP}{unemployment rate in ZIP code.}
#' }
#'
#' @source \url{http://people.stern.nyu.edu/wgreene/Microeconometrics.htm}
#'
#' @references
#' Greene, W. H. (1998).  Sample selection in credit-scoring models.
#' \emph{Japan and the World Economy}, 10(3), 299-316.
#'
#' @examples
#' \dontrun{
#' data(AmEx)
#' dim(AmEx)
#'
#' outcome   <- DEFAULT  ~ AGE + ACADMOS + ADEPCNT + AEMPMOS + MAJORDRG +
#'              MINORDRG + OWNRENT + INCOME + SELFEMPL
#' selection <- CARDHLDR ~ AGE + ACADMOS + ADEPCNT + AEMPMOS + MAJORDRG +
#'              MINORDRG + OWNRENT + INCOME + SELFEMPL + BANKSAV + BANKCH
#'
#' fit <- HeckSelect(selection, outcome, data = AmEx,
#'                    penalty = "alasso", Model = "AMH", crit = "bic")
#' summary(fit)
#' }
#'
#' @keywords datasets
"AmEx"