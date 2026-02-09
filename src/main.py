from fastapi import FastAPI

app = FastAPI(title="Biological Molecular Interaction Intelligence", version="0.1.0")

@app.get("/")
async def root():
    return {"message": "Welcome to the Biological Molecular Interaction Intelligence API"}

@app.get("/health")
async def health():
    return {"status": "ok"}
