#' Compare Two UMAP Images Using GPT-4o
#'
#' This function sends two UMAP plots to the OpenAI GPT-4o model and asks which one is visually more accurate,
#' based on cluster separation and boundary clarity. The function uses the OpenAI API.
#'
#' @param img1_path Path to the first image (PNG format).
#' @param img2_path Path to the second image (PNG format).
#' @param api_key Your OpenAI API key. If not provided, the function will use the environment variable \code{OPENAI_API_KEY}.
#' @param verbose Logical. If \code{TRUE}, prints the raw response from GPT-4o. Default is \code{FALSE}.
#'
#'
#' @return A named list:
#' \describe{
#'   \item{result}{An integer, 1 or 2, indicating the better UMAP image.}
#'   \item{raw_reply}{Character. The raw text reply from GPT-4o, useful for logging and debugging.}
#' }
#'
#' @import httr jsonlite base64enc
#'
#' @examples
#' \dontrun{
#' Sys.setenv(OPENAI_API_KEY = "sk-...")
#' result <- gpt_compare_umap("path/to/umap1.png", "path/to/umap2.png")
#' print(result$result)
#' }
#'
#' @export


gpt_compare_umap <- function(img1_path, img2_path, api_key = Sys.getenv("OPENAI_API_KEY"), verbose = FALSE) {

  if (is.null(api_key) || api_key == "" || startsWith(api_key, "your_openai_API_key")) {
    stop("OpenAI API key not found or invalid. Set via Sys.setenv(OPENAI_API_KEY = 'sk-...') or pass directly.")
  }

  if (!file.exists(img1_path) || !file.exists(img2_path)) {
    stop(paste("Image paths invalid:", img1_path, "or", img2_path))

  }

  img1_base64 <- base64enc::base64encode(img1_path)
  img2_base64 <- base64enc::base64encode(img2_path)
  img1_url <- paste0("data:image/png;base64,", img1_base64)
  img2_url <- paste0("data:image/png;base64,", img2_base64)

  message_body <- list(
    model = "gpt-4o",
    messages = list(
      list(
        role = "user",
        content = list(
          list(type = "text", text = paste(
            "You will be shown two UMAP plots.",
            "Please carefully examine both. Image 1 is labeled `1', and Image 2 is labeled `2'.",
            "Based on visual accuracy (cluster separation, boundary clarity, etc.), which one is better?",
            "Please respond with only: 1 or 2."

          )),
          list(type = "image_url", image_url = list(url = img1_url)),
          list(type = "image_url", image_url = list(url = img2_url))
        )
      )
    )
  )

  res <- httr::POST(
    url = "https://api.openai.com/v1/chat/completions",
    httr::add_headers(Authorization = paste("Bearer", api_key), `Content-Type` = "application/json"),
    body = jsonlite::toJSON(message_body, auto_unbox = TRUE)
  )

  result <- httr::content(res, as = "parsed", encoding = "UTF-8")

  if (!is.null(result$error)) {
    stop("OpenAI API Error: ", result$error$message)
  }
  if (is.null(result$choices[[1]]$message$content)) {
    stop("Invalid API response: GPT did not return a message.")
  }

  reply <- result$choices[[1]]$message$content
  if (verbose) message("GPT-4o raw reply: ", reply)

  decision <- if (grepl("\\b1\\b", reply)) 1 else if (grepl("\\b2\\b", reply)) 2 else NA

  if (is.na(decision)) {
    stop("GPT-4o returned unclear reply: ", reply)
  }

  return(list(result = decision, raw_reply = reply))
}

