#' Compare Two UMAP Images Using GPT-4o
#'
#' This function sends two UMAP plots to the OpenAI GPT-4o model and asks which one is visually more accurate,
#' based on cluster separation and boundary clarity. The function uses the OpenAI API.
#'
#' @param img1_path Path to the first image (PNG format).
#' @param img2_path Path to the second image (PNG format).
#' @param api_key Optional. Your OpenAI API key. If not provided, the function will use the environment variable `OPENAI_API_KEY`.
#'
#' @return A list with the selected image result (1 or 2) and the full GPT reply.
#' @export
#'
#' @import httr jsonlite base64enc
#'
#' @examples
#' \dontrun{
#' Sys.setenv(OPENAI_API_KEY = "sk-...")
#' result <- compare_umap("path/to/1.png", "path/to/2.png")
#' print(result$result)
#' print(result$raw_reply)
#' }
gpt_compare_umap <- function(img1_path, img2_path, api_key = Sys.getenv("OPENAI_API_KEY")) {
  if (api_key == "") {
    stop("OpenAI API key not found. Please set the `OPENAI_API_KEY` environment variable or pass `api_key` directly.")
  }

  if (!file.exists(img1_path) || !file.exists(img2_path)) {
    stop("One or both image paths are invalid or do not exist.")
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
          list(
            type = "text",
            text = paste(
              "You will be shown two UMAP plots. Please carefully examine both.",
              "Image 1 is labeled '1', and Image 2 is labeled '2'.",
              "Based on visual accuracy (cluster separation, boundary clarity, etc.), which one is better?",
              "Please respond with only: 1 or 2."
            )
          ),
          list(type = "image_url", image_url = list(url = img1_url)),
          list(type = "image_url", image_url = list(url = img2_url))
        )
      )
    )
  )

  attempt <- 1
  max_attempts <- 3
  success <- FALSE
  result <- NULL

  while (attempt <= max_attempts && !success) {
    res <- tryCatch({
      httr::POST(
        url = "https://api.openai.com/v1/chat/completions",
        httr::add_headers(
          Authorization = paste("Bearer", api_key),
          `Content-Type` = "application/json"
        ),
        body = jsonlite::toJSON(message_body, auto_unbox = TRUE)
      )
    }, error = function(e) e)

    if (inherits(res, "error")) {
      message(sprintf("Attempt %d failed: %s", attempt, res$message))
      attempt <- attempt + 1
      Sys.sleep(2)
    } else {
      result <- httr::content(res, as = "parsed", encoding = "UTF-8")
      if (!is.null(result$error)) {
        message(sprintf("API Error: %s", result$error$message))
        attempt <- attempt + 1
        Sys.sleep(2)
      } else {
        success <- TRUE
      }
    }
  }

  if (!success) {
    stop("Failed to get a valid response from OpenAI API after multiple attempts.")
  }

  reply <- result$choices[[1]]$message$content
  decision <- if (grepl("\\b1\\b", reply)) 1 else if (grepl("\\b2\\b", reply)) 2 else NA

  return(list(result = decision, raw_reply = reply))
}
